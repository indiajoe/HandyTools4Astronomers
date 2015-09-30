#!/usr/bin/env python
""" This is to fit Spectral Line models to reduced 1D spectrum. """

import sys
import os
import shlex
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.modeling import models, fitting

plt.rcParams['keymap.fullscreen']=[u'ctrl+f']  # To remove full screen togle on f keystrock

#LineToFit = 6563 #5577 #6300   # In wavelngth
FWindowToFit = 25  # In wavelgth, Forward size of window to fit the line model
BWindowToFit = 35  # In wavelgth, Forward size of window to fit the line model

SpectrumFileList = sys.argv[1]
OutputTableFile = sys.argv[2]
LinesToFitFile = sys.argv[3]


def LoadLinesToFit(filename):
    """ Loads line list form the file to fit the input spectrum 
    Format: 
    wavelength "Name" [Flag]

    where Flag is optional if 
         PCygni: fit double gaussian.
         BandPass:  give median value near bandpass value
    """
    OutLineList = []
    OutPCygniList = []
    OutBandPassList = []
    OutLineNameDic = dict()
    OutBandPassDic = dict()
    with open(filename) as linelistFile:
        linelist = [line.rstrip() for line in linelistFile if line[0] != '#']
    for line in linelist:
        Flag = 'Line'
        linesplit = shlex.split(line)
        if len(linesplit) > 2:
            Flag = linesplit[2]
        if Flag == 'Line':
            OutLineList.append(float(linesplit[0]))
            OutLineNameDic[float(linesplit[0])] = linesplit[1]
        elif Flag == 'PCygni':
            OutPCygniList.append(float(linesplit[0]))
            OutLineNameDic[float(linesplit[0])] = linesplit[1]
        elif Flag == 'BandPass':
            OutBandPassList.append(linesplit[1])
            OutBandPassDic[linesplit[1]] = float(linesplit[0])
            
    return OutLineList, OutPCygniList, OutBandPassList, OutLineNameDic, OutBandPassDic

def LoadSpectrumFile(filename,hdu=0,indx=0):
    """ Loads a returns the spectrum as 2 column numpy array """
    if os.path.splitext(filename)[-1] == '.npy':
        Spectrum = np.load(filename)
    elif os.path.splitext(filename)[-1] == '.fits':
        Spectrum = LoadFitsSpectrum(filename,hdu=hdu,indx=indx)
    else:
        print('Unrecognised input format')
        raise NotImplementedError()
    return Spectrum

def LoadFitsSpectrum(filename,hdu=0,indx=0):
    """Load and return the wavelength calibrated input HCT fits spectrum as 2 column numpy array. 
    hdu  : specifies the hdulist to read data and header
    indx  : specifies the column in the data to choose. In HCT, 0 for flux data and 2 for sky """
    fitsfile = fits.open(filename)
    flux = fitsfile[hdu].data  #[indx,0,:]
    w = WCS(fitsfile[hdu].header)
    Size = fitsfile[hdu].header['NAXIS1']
    
    # try :
    #     ref_pixel = fitsfile[hdu].header['CRPIX1']
    #     coord_ref_pixel = fitsfile[hdu].header['CRVAL1']
    #     wave_per_pixel = fitsfile[hdu].header['CDELT1']
    # except KeyError as e :
    #     print('Error: Missing keywords in fits header to do wavelength calibration')
    #     print(e)
    #     print('You might have entered wrong file name. Hence I am raising IOError')
    #     print('Enter the fits file name which is wavelength calibrated.')
    #     raise IOError
    # else:
    #    w_start=coord_ref_pixel - ((ref_pixel-1) * wave_per_pixel)  #Starting wavelength
    #    Wavelengths = w_start+np.arange(len(flux))*wave_per_pixel    


    CoordArray = np.zeros((Size,w.naxis))
    CoordArray[:,0] = np.arange(Size)
    Wavelengths = w.wcs_pix2world(CoordArray,0)[:,0]

    return np.vstack((Wavelengths,flux)).T


def NearestIndex(Array,value):
    """ Returns the index of element in numpy 1d Array nearest to value """
    return np.abs(Array-value).argmin()


def FitModel(X,Y,Model):
    """ Fits and returns the input astropy.model on X, Y data. """
    #Define fitting object
    fitter = fitting.LevMarLSQFitter()#SLSQPLSQFitter()
    #Fit the model to data
    Model_fit = fitter(Model, X, Y)  
    print(Model_fit)
    return Model_fit
    

def FitInteractiveGaussian(Spec,Line):
    """ Fits Interactively a gaussian, by asking user to select region on both sides to fit continuum """
    Start = NearestIndex(Spec[:,0], Line - BWindowToFit)
    End = NearestIndex(Spec[:,0], Line + FWindowToFit)

    #Initial Estimates
    Amp = np.max(np.abs(Spec[Start:End,1] - np.median(Spec[Start:End,1])))
    Bkgleft = [Start, Start+5]
    Bkgright = [End-5, End]

    #Define the line model to fit 
    LineModel = models.Gaussian1D(amplitude=Amp, mean=0, stddev=4)

    # Background straight line polynomial
    Bkg = np.poly1d(np.polyfit([np.mean(Spec[Bkgleft[0]:Bkgleft[1],0]),np.mean(Spec[Bkgright[0]:Bkgright[1],0])],\
                           [np.median(Spec[Bkgleft[0]:Bkgleft[1],1]),np.median(Spec[Bkgright[0]:Bkgright[1],1])], 1))

    # Fit the model Only inside the region excluding the Bkg estimation region
    BkgEstimate = Bkg(Spec[Bkgleft[1]:Bkgright[0],0])
    FittedModel = FitModel(Spec[Bkgleft[1]:Bkgright[0],0]-Line, Spec[Bkgleft[1]:Bkgright[0],1]-BkgEstimate,LineModel)  # Centering the x on zero for fitting

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('q<->w to select left Bkg.  e<->r to select right Bkg. f to (re)fit.')
    ax.plot(Spec[Start:End,0]-Line,Spec[Start:End,1],linestyle='--', drawstyle='steps-mid',marker='.')
    bkgLine, = ax.plot(Spec[Bkgleft[1]:Bkgright[0],0]-Line,BkgEstimate,color='g')
    GaussLine, = ax.plot(Spec[Bkgleft[1]:Bkgright[0],0]-Line,BkgEstimate+FittedModel(Spectrum[Bkgleft[1]:Bkgright[0],0]-Line),color='r')
    PlotedLines = [bkgLine,GaussLine]
    ModelBkgBkgList = [FittedModel,BkgEstimate,Bkg]
    # Define the function to run while key is pressed
    def on_key(event):
        if event.key == 'q' :
            Bkgleft[0] = NearestIndex(Spec[:,0], Line + event.xdata)
            print 'Left Continuum L Position = {0}, Index = {1}'.format(event.xdata, Bkgleft[0] )
        elif event.key == 'w' :
            Bkgleft[1] = NearestIndex(Spec[:,0], Line + event.xdata)
            print 'Left Continuum R Position = {0}, Index = {1}'.format(event.xdata, Bkgleft[1] )
        elif event.key == 'e' :
            Bkgright[0] = NearestIndex(Spec[:,0], Line + event.xdata)
            print 'Right Continuum L Position = {0}, Index = {1}'.format(event.xdata, Bkgright[0] )
        elif event.key == 'r' :
            Bkgright[1] = NearestIndex(Spec[:,0], Line + event.xdata)
            print 'Right Continuum R Position = {0}, Index = {1}'.format(event.xdata, Bkgright[1] )

        elif event.key == 'f' :
            #Sanity Check
            if (len(Spec[Bkgleft[0]:Bkgleft[1],1]) > 0) and (len(Spec[Bkgright[0]:Bkgright[1],1]) > 0) and (Bkgright[0] > Bkgleft[1]):
                
                Bkg = np.poly1d(np.polyfit([np.mean(Spec[Bkgleft[0]:Bkgleft[1],0]),np.mean(Spec[Bkgright[0]:Bkgright[1],0])],\
                                 [np.median(Spec[Bkgleft[0]:Bkgleft[1],1]),np.median(Spec[Bkgright[0]:Bkgright[1],1])], 1))

                BkgEstimate = Bkg(Spec[Bkgleft[1]:Bkgright[0],0])
                FittedModel = FitModel(Spec[Bkgleft[1]:Bkgright[0],0]-Line, Spec[Bkgleft[1]:Bkgright[0],1]-BkgEstimate,LineModel)
                ModelBkgBkgList[0],ModelBkgBkgList[1],ModelBkgBkgList[2] = FittedModel,BkgEstimate,Bkg

                ax.lines.remove(PlotedLines[0])  #Removing previous line plots
                ax.lines.remove(PlotedLines[1])  

                PlotedLines[0], = ax.plot(Spec[Bkgleft[1]:Bkgright[0],0]-Line,BkgEstimate,color='g')
                PlotedLines[1], = ax.plot(Spec[Bkgleft[1]:Bkgright[0],0]-Line,BkgEstimate+FittedModel(Spectrum[Bkgleft[1]:Bkgright[0],0]-Line),color='r')
                ax.figure.canvas.draw()
            else:
                print('Error: Incompatible Background estimation positions given by user')

    cid = fig.canvas.mpl_connect('key_press_event', on_key)
    plt.show()

    FittedModel,BkgEstimate,Bkg = ModelBkgBkgList[0],ModelBkgBkgList[1],ModelBkgBkgList[2]
    
    Area = FittedModel.amplitude.value * abs(FittedModel.stddev.value) *np.sqrt(np.pi)
    Eqw = Area/Bkg(FittedModel.mean.value+Line)
    LineBkgSub = Spec[Bkgleft[1]:Bkgright[0],1] - Bkg(Spec[Bkgleft[1]:Bkgright[0],0])
    DeltaX = np.abs(np.mean(Spec[Bkgleft[1]:Bkgright[0]-1,0] - Spec[Bkgleft[1]+1:Bkgright[0],0]))
    SumAbove = np.sum(LineBkgSub[LineBkgSub>0]) * DeltaX
    SumBelow = np.sum(LineBkgSub[LineBkgSub<0]) * DeltaX
    print('Flux = {0} ; EQW ={1}'.format(Area/FLUXSCALE,Eqw*-1))
    print('Flux Above= {0} ; Flux Below ={1}'.format(SumAbove/FLUXSCALE,SumBelow/FLUXSCALE))
    return FittedModel,BkgEstimate,Area,Eqw,SumAbove,SumBelow

def FitSingleGaussianLine(Spec,Line):
    """ Fits a stright line background model and Single Gaussian """
    Start = NearestIndex(Spec[:,0], Line - BWindowToFit)
    End = NearestIndex(Spec[:,0], Line + FWindowToFit)

    #Initial Estimates
    BackG = np.median(Spec[Start:End,1])    
    Amp = np.max(np.abs(Spec[Start:End,1] - BackG))

    #Define the line model to fit 
    LineModel = models.Gaussian1D(amplitude=Amp, mean=0, stddev=4) + models.Linear1D(slope=0,intercept=BackG)

    FittedModel = FitModel(Spec[Start:End,0]-Line, Spec[Start:End,1],LineModel)  # Centering the x on zero for fitting
    Area = FittedModel.amplitude_0.value * abs(FittedModel.stddev_0.value) *np.sqrt(np.pi)
    Eqw = Area/FittedModel[1](FittedModel.mean_0.value)
    print('Flux = {0} ; EQW ={1}'.format(Area/FLUXSCALE,Eqw))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(Spec[Start:End,0], Spec[Start:End,1]/FLUXSCALE, 'ko')
    ax1.plot(Spec[Start:End,0], FittedModel(Spectrum[Start:End,0]-Line)/FLUXSCALE, 'r-', lw=2)
    ax1.axvline(x=Line,linestyle='--',color='k')
    ax1.set_xlabel(r'$\lambda$ $(\dot{A})$')
    ax1.set_ylabel('$F_\lambda$')
    c = 2.998e5  # Light speed in km/s
    VelocityTicks = np.arange(-2000,1500,200)
    VelocityTicksLocation = VelocityTicks*Line/c + Line
    ax2.set_xticks(VelocityTicksLocation)
    ax2.set_xticklabels(VelocityTicks)
    ax2.set_xbound(ax1.get_xbound())
#    ax2.plot(c*(Spec[Start:End,0]-Line)/Line, Spectrum[Start:End,1]/FLUXSCALE,alpha=0) # Create a dummy plot
    ax2.grid(True)
    ax2.set_xlabel('$v$ $(km/s)$')
    plt.show(block=False)
    FitQuality = raw_input('Enter quality of fit (0=Good,1=Poor,2=Bad,3=Wrong etc..) :').strip()
    plt.close() 
    return FittedModel, FitQuality

def FitPCygniLine(Spec,Line,interactive=False):
    """ Fits a straight line background model and two Single Gaussians """
    Start = NearestIndex(Spec[:,0], Line - BWindowToFit)
    End = NearestIndex(Spec[:,0], Line + FWindowToFit)

    #Initial Estimates
    BackG = np.median(Spec[Start:End,1])    
    EmiAmp = [np.max(np.abs(Spec[Start:End,1] - BackG))]
    AbsAmp = [-1*EmiAmp[0]/2.0]
    EmiPos = [0]
    AbsPos = [-9]


    if interactive:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Left click Absoption and Right click Emission peaks')
        ax.plot(Spec[Start:End,0]-Line,Spec[Start:End,1])
        # Define the function to run while clicked.
        def onclick(event):
            if event.button == 1 :
                AbsPos[0], AbsAmp[0] = event.xdata, event.ydata - BackG
                print 'Absorb Position = {0}, Amplitude = {1}'.format(AbsPos[0], AbsAmp[0])
            elif event.button == 3 :
                EmiPos[0], EmiAmp[0] = event.xdata, event.ydata - BackG
                print 'Emission Position = {0}, Amplitude = {1}'.format(EmiPos[0], EmiAmp[0])

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
                

    #Define the line model to fit 
    LineModel = models.Gaussian1D(amplitude=EmiAmp[0], mean=EmiPos[0], stddev=3) \
                        + models.Gaussian1D(amplitude=AbsAmp[0], mean=AbsPos[0], stddev=3) \
                        + models.Linear1D(slope=0,intercept=BackG)

    # Limiting the gaussian stddev to be less than 5
    LineModel.stddev_0.min = -5
    LineModel.stddev_1.min = -5
    LineModel.stddev_0.max = 8
    LineModel.stddev_1.max = 8
    LineModel.amplitude_0.min = 0
    LineModel.amplitude_1.max = 0
    

    FittedModel = FitModel(Spec[Start:End,0]-Line, Spec[Start:End,1],LineModel)  # Centering the x on zero for fitting

    Area0 = FittedModel.amplitude_0.value * abs(FittedModel.stddev_0.value) *np.sqrt(np.pi)
    Eqw0 = Area0/FittedModel[2](0.0)
    Area1 = FittedModel.amplitude_1.value * abs(FittedModel.stddev_1.value) *np.sqrt(np.pi)
    Eqw1 = Area1/FittedModel[2](0.0)

    print('Flux_0 = {0} ; EQW_0 ={1}'.format(Area0/FLUXSCALE,Eqw0))
    print('Flux_1 = {0} ; EQW_1 ={1}'.format(Area1/FLUXSCALE,Eqw1))

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(Spec[Start:End,0], Spec[Start:End,1]/FLUXSCALE, 'ko')
    ax1.plot(Spec[Start:End,0], FittedModel(Spectrum[Start:End,0]-Line)/FLUXSCALE, 'r-', lw=2)
    ax1.plot(Spec[Start:End,0], (FittedModel[0](Spectrum[Start:End,0]-Line)+FittedModel[2](Spectrum[Start:End,0]-Line))/FLUXSCALE, 'g-', lw=1.5)
    ax1.plot(Spec[Start:End,0], (FittedModel[1](Spectrum[Start:End,0]-Line)+FittedModel[2](Spectrum[Start:End,0]-Line))/FLUXSCALE, 'g-', lw=1.5)
    ax1.axvline(x=Line,linestyle='--',color='k')
    ax1.set_xlabel(r'$\lambda$ $(\dot{A})$')
    ax1.set_ylabel('$F_\lambda$')
    c = 2.998e5  # Light speed in km/s
    VelocityTicks = np.arange(-2000,1500,200)
    VelocityTicksLocation = VelocityTicks*Line/c + Line
    ax2.set_xticks(VelocityTicksLocation)
    ax2.set_xticklabels(VelocityTicks)
    ax2.set_xbound(ax1.get_xbound())
    ax2.grid(True)
    ax2.set_xlabel('$v$ $(km/s)$')
    plt.show(block=False)
    FitQuality = raw_input('Enter quality of fit (0=Good,1=Poor,2=Bad,3=Wrong etc..) :').strip()
    plt.close() 
    return FittedModel, FitQuality



##################################################
LinesToFitGaussian, LinesToFitPCygni, BandPassFiltersToMeasure, LineNameDic, FiltWL = LoadLinesToFit(LinesToFitFile)
# BandPassFiltersToMeasure = [] #['J','H','Ks']
# LinesToFitPCygni = []
# LinesToFitGaussian = [12818] # [10829,10938]#, JOS

DATEOBSHDR = 'DATE'  # Fits header for Date

# Dictinary of line wavelengths and element for later reference
# LineNameDic={3889.0:'H8',3970.1:'H7',4101.7:'Hd',4340.5:'Hg',4861.3:'Hb',3933.6:'CaII(K)',4069:'[SII]',4815:'[FeII]',5015.6:'HeI',5895:'NaD',6300.3:'[OI]',6363.8:'[OI]',6432:'FeII',6517:'FeII',6562.8:'Ha',6730:'[SiII]',7155:'[FeII]',7291:'[CaII]',7324:'[OI]+[CaII]',7380:'[FeII]',7699:'KI',7773:'OI',8388:'FeI',8446:'OI',8616:'[FeII]',8498:'CaII',8542:'CaII',8662:'CaII',5577:'[OI]',10829:'HeI',10938:'PaG',12818:'PaB'}

# FiltWL = {'V':5500,'R':6400,'I':7900} # Johnson V, Cousins R,I  J H Ks etc.. if BandPassFiltersToMeasure has it...



with open(SpectrumFileList) as listfile:
    SpectrumList=[fls.rstrip() for fls in listfile if fls[0] != '#']


ListOfTableRows=[]
for SpectrumFile in SpectrumList:
    print('Line Fitting on File : {0}'.format(SpectrumFile))
    # Initialise dictionary for the file
    TableRow = {'File':SpectrumFile}#,'Date':fits.getval(SpectrumFile,DATEOBSHDR)}#,'JD':fits.getval(SpectrumFile,'JD')}  
    TableColumns = ['File']#,'Date']#,'JD']
    # # First fit sky lines [OI] 5577.3 and [OI] 6300.3 to find wavelength shift correction    
    # ## Load the sky spectrum first
    # Spectrum = LoadFitsSpectrum(SpectrumFile,indx=2)
    # # Scale the spectrum to good number
    # FLUXSCALE = 1.0/np.median(Spectrum[:,1]) #10**15
    # Spectrum[:,1] *= FLUXSCALE
    # print('Raw Median Flux ={0}'.format(FLUXSCALE))

    # Wshifts=[]
    # for line in [5577.3,6300.3]:
    #     FittedLine, Quality = FitSingleGaussianLine(Spectrum,line)
    #     if Quality in ['0','1'] :Wshifts.append(FittedLine.mean_0.value)

    # WL_SHIFT_CORR = np.mean(Wshifts)
    # print('Wavelength correction shifts: {0} \n Mean correction = {1}'.format(Wshifts,WL_SHIFT_CORR))
    # print('STDEV of the shifts calulated = {0}'.format(np.std(Wshifts)))

    FLUXSCALE = 1
    WL_SHIFT_CORR = 0
    ## Load the Flux spectrum now
    Spectrum = LoadSpectrumFile(SpectrumFile,hdu=0,indx=0)
    # Scale the spectrum to good number and also correct the wavelength shift
    Spectrum[:,1] *= FLUXSCALE
    Spectrum[:,0] -= WL_SHIFT_CORR

    TableRow['WLshiftCorr'] = 0#WL_SHIFT_CORR
    TableColumns.append('WLshiftCorr')

    ## V,R,I band fluxes from spectrum..

    for filt in BandPassFiltersToMeasure : #['V','R','I']: #]: # For G8: 
        FluxM =  np.median(Spectrum[NearestIndex(Spectrum[:,0], FiltWL[filt]-BWindowToFit): NearestIndex(Spectrum[:,0], FiltWL[filt] + FWindowToFit),1])/FLUXSCALE
        print('Band :{0}  F_lambda = {1}'.format(filt,FluxM))
        TableRow['Flux_'+filt] = FluxM
        TableColumns.append('Flux_'+filt)

    # Now start fitting lines in spectrum
    print('#>>> Start Fitting PCygni profiles....')
    for line in LinesToFitPCygni : #[6300.3,6562.8]: #[4861.3,6300.3,6562.8]:#  For G8 :[6300.3,6562.8,8542,8662]
        print('*'*6+'{0}'.format(line)+'*'*6)
        FittedLine, Quality = FitPCygniLine(Spectrum,line,interactive=True)
        if Quality in ['0','1','2']:
            AreaEmi = FittedLine.amplitude_0.value * abs(FittedLine.stddev_0.value) *np.sqrt(np.pi)
            EqwEmi = AreaEmi/FittedLine[2](FittedLine.mean_0.value)
            AreaAbs = FittedLine.amplitude_1.value * abs(FittedLine.stddev_1.value) *np.sqrt(np.pi)
            EqwAbs = AreaAbs/FittedLine[2](FittedLine.mean_1.value)
            PeakPosEmi = FittedLine.mean_0.value + line
            PeakPosAbs = FittedLine.mean_1.value + line
        else:
            print('*'*5+'Fitting Peaks and Dips with seperate gaussians, close window to continue')
            print('>>>PLEASE: First fit the ABSORPTION component alone')

            FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
            AreaAbs = Area
            EqwAbs = Eqw 
            PeakPosAbs = FittedLine.mean.value + line

            print('>>>PLEASE: Now fit the EMISSION component alone')
            FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
            AreaEmi = Area
            EqwEmi = Eqw 
            PeakPosEmi = FittedLine.mean.value + line
            

        TableRow['FluxGPEmi_'+str(line)] = AreaEmi/FLUXSCALE
        TableRow['FluxGPAbs_'+str(line)] = AreaAbs/FLUXSCALE
        TableRow['eqwGPEmi_'+str(line)] = EqwEmi *-1
        TableRow['eqwGPAbs_'+str(line)] = EqwAbs *-1
        TableRow['PeakGPEmi_'+str(line)] = PeakPosEmi
        TableRow['PeakGPAbs_'+str(line)] = PeakPosAbs
        TableRow['FitQPCyg_'+str(line)] = Quality
        TableColumns+=['FluxGPEmi_'+str(line),'FluxGPAbs_'+str(line),'eqwGPEmi_'+str(line),'eqwGPAbs_'+str(line),'PeakGPEmi_'+str(line),'PeakGPAbs_'+str(line),'FitQPCyg_'+str(line)]
 
        ### Uncomment the lines below to do manual fit seperately each time....
        # print('*'*5+'Fitting Peaks and Dips with seperate gaussians, close window to continue')
        # print('>>>PLEASE: First fit the Absorbtion component alone')
        # FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
        # TableRow['FluxMGPAbs_'+str(line)] = Area/FLUXSCALE 
        # TableRow['eqwMGPAbs_'+str(line)] = Eqw *-1
        # TableRow['PeakMGPAbs_'+str(line)] = FittedLine.mean.value + line
        # TableColumns+=['FluxMGPAbs_'+str(line),'eqwMGPAbs_'+str(line),'PeakMGPAbs_'+str(line)]
        # print('>>>PLEASE: Now fit the Emission component alone')
        # FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
        # TableRow['FluxMGPEmi_'+str(line)] = Area/FLUXSCALE 
        # TableRow['eqwMGPEmi_'+str(line)] = Eqw *-1
        # TableRow['PeakMGPEmi_'+str(line)] = FittedLine.mean.value + line
        # TableColumns+=['FluxMGPEmi_'+str(line),'eqwMGPEmi_'+str(line),'PeakMGPEmi_'+str(line)]

        print('*'*15)

    print('#>>> Start Fitting Single Gaussian profiles....')
    for line in  LinesToFitGaussian : #[6300.3,6363.8,6432,6517,6562.8,6730,7155,7291,7324,7380,7699,7773]:# [3889.0,3970.1,4101.7,4340.5,4861.3,3933.6,4069,4815,5015.6,5895,6300.3,6363.8,6432,6517,6562.8,6730,7155,7291,7324,7380,7699,5577]:#For G8 :[5895,6300.3,6363.8,6432,6517,6562.8,6730,7155,7291,7324,7380,7699,7773,8388,8446,8616,8498,8542,8662,5577]
        print('*'*6+'{0}'.format(line)+'*'*6)
        FittedLine, Quality = FitSingleGaussianLine(Spectrum,line)
        if Quality in ['0','1','2']:
            Area = FittedLine.amplitude_0.value * abs(FittedLine.stddev_0.value) *np.sqrt(np.pi)
            Eqw = Area/FittedLine[1](FittedLine.mean_0.value)
            PeakPos = FittedLine.mean_0.value + line
        else:
            print('*'*5+'Fitting gaussian manually, close window to continue')
            FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
            PeakPos = FittedLine.mean.value + line
            
        TableRow['FluxG_'+str(line)] = Area/FLUXSCALE
        TableRow['eqwG_'+str(line)] = Eqw *-1
        TableRow['PeakG_'+str(line)] = PeakPos
        TableRow['FitQ_'+str(line)] = Quality
        TableColumns+=['FluxG_'+str(line),'eqwG_'+str(line),'PeakG_'+str(line),'FitQ_'+str(line)]

        ### Uncomment the lines below to do manual fit seperately each time....
        # print('*'*5+'Fitting gaussians manually, close window to continue')
        # FittedLine,BkgEstimate,Area,Eqw,SumAbove,SumBelow = FitInteractiveGaussian(Spectrum,line)
        # TableRow['FluxMG_'+str(line)] = Area/FLUXSCALE 
        # TableRow['eqwMG_'+str(line)] = Eqw *-1
        # TableRow['PeakMG_'+str(line)] = FittedLine.mean.value + line
        # TableRow['FluxMSumEmi_'+str(line)] = SumAbove/FLUXSCALE 
        # TableRow['FluxMSumAbs_'+str(line)] = SumBelow/FLUXSCALE 
        # TableColumns+=['FluxMG_'+str(line),'eqwMG_'+str(line),'PeakMG_'+str(line),'FluxMSumEmi_'+str(line),'FluxMSumAbs_'+str(line)]

        print('*'*15)
    
    # Append the this table row
    ListOfTableRows.append(TableRow)

# Generate an astropy table
OutputTable = Table(rows=ListOfTableRows,names=TableColumns)  # Giving names to keep the column orders correct
OutputTable.show_in_browser(jsviewer=True) 
OutputTable.write(OutputTableFile,format='ascii.fixed_width')


