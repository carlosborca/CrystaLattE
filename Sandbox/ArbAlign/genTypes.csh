#!/bin/tcsh

if ( $#argv != 1 ) then
   echo "Error in the number of arguments passed"
   echo "Usage: genTypes.csh <XYZ file name>"
   echo "   Eg: genTypes.csh test.xyz"
   exit
endif

#source /usr/local/Modules/3.2.10/init/tcsh
#module load openbabel/default

head -2 $1 > t-$1
babel -ixyz $1 -osy2 |grep -A 10000 "@<TRIPOS>ATOM" |awk '{if($1=="@<TRIPOS>BOND") {exit} else if(NF>3) {print $6, $3, $4, $5} else {} }' >> t-$1 
