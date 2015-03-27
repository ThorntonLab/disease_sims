## ----,fig.width=5,fig.height=5-------------------------------------------
library(diseaseSims)
idx=read.table(system.file("extdata", "index.txt", package = "diseaseSims"))
popfile=system.file("extdata", "pop.bin.gz", package = "diseaseSims")
OFFSET=0
for(i in 1:nrow(idx))
    {
      ##Get the risk variant matrix for the population
      XX=getRiskVariantMatrix(popfile,OFFSET,model="additive")
      XX.aov = vpv1aov(XX)
      if( i == 1 )
      {
	plot(XX.aov[,1],XX.aov[,2],
		xlab="Mutation frequency",
		ylab="Variance explained",
		xlim=c(0,0.1),type="l",
		main="3 replicates, additive model")
      }
      else
      {
	lines(XX.aov[,1],XX.aov[,2],col=i)
      }
      OFFSET=OFFSET+idx$V4[i]
    }

