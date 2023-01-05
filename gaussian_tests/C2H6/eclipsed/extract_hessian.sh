#/bin/bash

fname=$(ls *.chk)
fbname=$(basename $fname .chk)

#echo "$fname $fbname" 

formchk $fname 
cat << EOF > demofc.inp 
$fbname.fchk
Atomic numbers
Current cartesian coordinates
Cartesian Force Constants
quit 
EOF
demofc < demofc.inp > hessian.out 

rm $fbname.fchk demofc.inp 

