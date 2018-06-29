#!/bin/bash


par=$1
tim=$2

psr=$(grep PSRJ $par | awk '{print $2}')

./disableparam.py $par START FINISH > $par.tmp
./enableparam.py $par.tmp JUMP RAJ DECJ F0 F1 PX PB A1 OM ECC T0 -glitch > $par.tmp2
cat $par.tmp2 | grep -v PM | grep -v F2 > $par.m
rm $par.tmp $par.tmp2

#tempo2 -f $par.m $tim -qrfit -nofit -fit pmra -fit pmdec -fit f2 -outpar $par.f2

#f2range=$(grep F2 $par.f2 | awk '{print sqrt($2*$2)+8*$4}')
#pmra_r=$(grep PMRA $par.f2 | awk '{print sqrt($2*$2)+8*$4}')
#pmdec_r=$(grep PMRA $par.f2 | awk '{print sqrt($2*$2)+8*$4}')

f0=$(grep -e '^F0 ' $par | awk '{print $2}')
f1=$(grep -e '^F1 ' $par | awk '{print $2}')

f2range=$(echo $f0 $f1 | awk '{print(25000*$2*$2/$1)}')


pmra_r=1000

echo "F2 range is: $f2range"

redmax=-9

tempo2 -f $par.m $tim -qrfit -outpar ${psr}_run.par

../pmtest/run_enterprise.py ${psr}_run.par $tim -j --all-corner --red-ncoeff 60 -N 1e5 -M "--Ared-max $redmax --pm --pm-range $pmra_r --f2 $f2range" \
    "--Ared-max $redmax --pm --pm-range $pmra_r" \
    "--no-red-noise --pm --pm-range $pmra_r --f2 $f2range" \
    "--no-red-noise --pm --pm-range $pmra_r"

