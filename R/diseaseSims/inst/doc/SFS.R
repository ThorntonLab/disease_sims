## ---- fig.height=5, fig.width=5------------------------------------------
library(diseaseSims)
library(dplyr)
##n = 100, and if you use n >= 2N, 2N will be used
##seed = 101
sfs.all = sfs(system.file("extdata", "pop.bin.gz", package = "diseaseSims"),100,101)

##Normalize the SFS
sfs.norm = sfs.all %>%
    group_by(replicate) %>%
        mutate( n.norm = n/sum(n),
                r.norm = r/sum(r) )

##Average, normalized SFS
sfs.norm.mean = sfs.norm %>%
    group_by(i) %>%
        summarise(n.norm.mean = mean(n.norm),
                  r.norm.mean = mean(r.norm))

plot(sfs.norm.mean$i,sfs.norm.mean$n.norm.mean,xlab=("Frequency (out of 100)"),
								ylab="Proportion of mutations",
								ylim=c(0,0.5),
								pch=17)
points(sfs.norm.mean$i,sfs.norm.mean$r.norm.mean,pch=18,col="red")
legend("topright",c("Neutral","Risk"),pch=c(17,18),col=c("black","red"),bty="n")

