#!sh

for i in 0.01 0.025 0.05 0.075 0.1 0.125 0.175 0.25 0.35 0.5
do
    cd lambda$i
    qsub -N SIM$i ../runtest.sh `pwd` $i
    qsub -N MAKECC$i -hold_jid SIM$i ../makeCC.sh `pwd`
    qsub -N SMARKER$i -hold_jid MAKECC$i ../singleMarker.sh `pwd`
    qsub -N ESM$i -hold_jid MAKECC$i ../esmTest.sh `pwd`
    cd ..
done