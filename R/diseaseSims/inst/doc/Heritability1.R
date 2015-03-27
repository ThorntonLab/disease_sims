## ------------------------------------------------------------------------
idx=read.table(system.file("extdata", "index.txt", package = "diseaseSims"))
##open in binary mode...
f=gzfile(system.file("extdata", "phenotypes.bin.gz", package = "diseaseSims"),"rb")

h2=array()
for(i in 1:nrow(idx))
{
	N = readBin(f,"integer",1)
	m = matrix(readBin(f,"numeric",2*N),ncol=2,byrow=TRUE)
	##column 1 = genetic component, column2 = random component
	h2[i]=var(m[,1])/( var(m[,1]) + var(m[,2]) )
}
##Now, you can get the mean, etc., which will not be very accurate for this example,
##as it is based on just 3 replicates
print(mean(h2))
print(h2)

