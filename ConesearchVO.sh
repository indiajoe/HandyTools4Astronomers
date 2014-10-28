#!/bin/bash
#---------------
# This is a simple script to download cone search results from a catalog using VO.
# Usage: ConesearchVO.sh InputFile.txt
#---------------
# Format of columns in InputFile.txt should be as follows
# RA  DEC  Radius  Catalog

# Where RA, DEC and Radius should be in degrees.
# Catalog is the VO mname of the catalog. Commonly used ones are given below.
# wise_allwise_p3as_psd  -> AllWISE Source Catalog
# fp_psc                 -> 2MASS Point Source Catalog
# glimpse_s07            -> GLIMPSE I Spring 07 Catalog (Spitzer)
# cosmos_phot            -> COSMOS Photometry Catalog
# iraspsc                -> IRAS Point Source Catalog 
#-------------------------------------Enjoy!! indiajoe

OUTformat='csv'  # Available output formats => votable, ipac_table, csv, tsv, fits

# Sanity check
if [ ! -f $1 ]; then
    echo "Input file $1 not found"
    echo "Usage: $0 InputFile.txt"
    exit 1
fi

while read -a line
do
    RA=${line[0]}
    DEC=${line[1]}
    Radius=${line[2]}
    Catalog=${line[3]}
    OutFile=${RA//./}${DEC//./}.$OUTformat
    echo "Quering $RA $DEC in radius $Radius from Catalog $Catalog"
    echo "Output will be saved to $OutFile"
    wget -O "$OutFile" "http://irsa.ipac.caltech.edu/SCS?table=$Catalog&RA=$RA&DEC=$DEC&SR=$Radius&format=$OUTformat" 
done < $1

# That is all. For more details on cone search via http by VO visit
# http://irsa.ipac.caltech.edu/docs/vo_scs.html
 
