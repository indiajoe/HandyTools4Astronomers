#!/bin/bash
#This shell script is to print the fits files with specific key values in header
if [ $# -lt 3 ] ; then
    echo
    echo " This script will list the fits files with specific key values in header"
    echo " Usage :"
    echo " $0 HEADERKEY VALUE *.fits "
    echo " It is important to have all the above entires in same order."
    echo " Enjoy !!!! ----------indiajoe"
    echo
    exit 2
fi
headerkey=$1
shift
value=$1
shift
for file in "$@"; do
    HeaderCard=$(fold $file | sed '/^END/q' | grep $headerkey)
    if [ $? -ne 0 ]; then
	echo "ERROR: Header keyword $headerkey missing in the fits file $file" 1>&2
    else 
	ValueComment=${HeaderCard#*=}  #Removing everything before =
	keyValue=${ValueComment%/*}   #Removing everyithing after /
	if [[ $keyValue =~ $value ]] ; then
	    echo $file
	fi
    fi
done

