library(Rcpp)
sourceCpp("Rcpp_helper.cc")
#f="/dfs1/test/bio/krthornt/does_model_matter/dist/nogrowth/additive/lambda0.05/pop.bin.gz"
#readPop(f,622932,12)
f="/dfs1/test/bio/krthornt/does_model_matter/dist/nogrowth/recessive/lambda0.5/effects.bin.gz"
e=getEsizes(f,767220)
f="/dfs1/test/bio/krthornt/does_model_matter/dist/nogrowth/recessive/lambda0.5/ccfile.bin"
cc=getCCblock(f,141481544)

#print(e[1:10,])
#print(cc$pos)
f="/dfs1/test/bio/krthornt/does_model_matter/dist/nogrowth/recessive/lambda0.5/phenotypes.bin.gz"
p=getPheno(f,37120464)

ccp=rowSums(cc$phenos)
pp=rowSums(p)

for(i in 1:length(ccp))
    {
        print(which(pp==ccp[i]))
    }
