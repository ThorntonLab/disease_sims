#!sh

#$ -q krt,bio

cd $1

module load R

echo "esize reseq gwas"
for i in 0.01 0.025 0.05 0.075 0.1 0.125 0.175 0.25 0.35 0.5
do
    z=`R --no-save --slave --args lambda$i/logit_pvalues_index.txt lambda$i/logit_pvalues.txt < ~/src/disease_sims/logit_power.R|cut -d" " -f2,3 | perl -p -i -e 's/\"//go'`
    echo $i $z
done