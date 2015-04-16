#!/bin/bash

parallel -k echo ::: {1..300} | parallel --jobs 4 R_LIBS=$HOME Rscript burdenPerm.R 0 ccindex.txt ccfile.bin.gz perms.gz 10000
R_LIBS=$HOME Rscript burdenPermP.R 0 ccindex.txt ccfile.bin.gz perms.gz pvals.txt.gz
