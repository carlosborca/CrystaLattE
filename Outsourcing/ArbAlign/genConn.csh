#!/bin/tcsh

if ( $#argv != 1 ) then
   echo "Error in the number of arguments passed"
   echo "Usage: genConn.csh <XYZ file name>"
   echo "   Eg: genConn.csh test.xyz"
   exit
endif

#source /usr/local/Modules/3.2.10/init/tcsh
#module load openbabel/default

babel -ixyz $1 -omna -xL1 | sed "s/[()-]//g" > .temp 
paste .temp $1 | awk '{if(NR<3){print $NF} else {printf "%s-%s\t %-f\t %-f\t %-f\n", $2,$1, $3, $4,$5}}' > c-$1
rm -f .temp

