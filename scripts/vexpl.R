library(diseaseSims)
##Silence the loading of dplyr...
suppressMessages(library(dplyr))

##Read in simulation index
idx=read.table("simindex.txt",colClasses=c("integer",rep("numeric",3)))

##Calculate the risk variant matrix under the recessive trait model (e.g., the same one we simulate under...) 
risk.matrix = getRiskVariantMatrix("pop.bin.gz",procIndex(idx,4,0),model="recessive")
##Fit the additive model using base::lm()
risk.matrix.aov = vpv1aov(risk.matrix)

pdf("vexpl.pdf",height=10,width=10,pointsize=18)
plot(risk.matrix.aov$p,risk.matrix.aov$rsq,
     type = "l",
     xlab="Risk mutation freq.",
     ylab="Cumulative % VG explained by additive effects",
     main="")
dev.off()
