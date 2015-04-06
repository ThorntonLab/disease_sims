# An example of all the steps

This folder contains an executable bash script illustrating a complete workflow involving a single simulated replicate.

To execute the script:

~~~{sh}
./example_workflow.sh
~~~

In order for this to work, you need ALL of the dependencies installed first!!! (Yes, there are a lot of them.)

This script probably takes 30-45 minutes to complete on a typical workstation.

## Notes

### Where to go for more info
The details of the R scripts show how various features of the diseaseSims R/Rcpp package work.  See the main documentation for that package for more detail.

### Translating these scripts to a cluster
In practice, we do not use workflows that are this simple.  Instead, we automate the processing of large batches of replicates via array jobs submitted to a cluster computing environment.

The main things to keep track of when using array jobs are:

* Keeping track of each replicate.  The command line programs typically take an argument to specify an integer that will serve as a unique ID for a replicate.  These ID numbers are recorded in the index file output by that program.  This scheme ensures that you can match up simulated population X with case/control panel X, etc.  The R package provides functions for processing these index files and reading in the correct record from the files.  
* Keeping track of the model.  In general, if you simulate data under an additive model, you probably want to analyze it under the same model.  There are, of course, exceptions to this rule.  For example, we apply logistic regressions under additive models to all simulated data, because that it how GWAS are done in practice.