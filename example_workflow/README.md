#How to redo an entire paper in a few commands

These scripts recreate the simulations from Thornton, Foran and Long (2013) PMID 23437004.

They are currently set up to run on the UCI cluster, which uses a Grid Engie (GE) queue.

To run these on your own GE system, you will need to:

* Edit "submit.sh" and add a -q queue1,queue2,etc. flag to each qsub command.
* Edit each script so that the path to each binary is correct

To execute the workflow:
```
sh setup.sh
sh submit.sh
```

When the simulations are done, the power of a single-marker GWAS (logistic regression of case/control status onto genotype under an additive model is obtained via

```
sh logit_power_table . > logit_power.txt
```

executed from the folder containing all the lambda\* subfolders.  (__NOTE:__ you may need to edit paths in that shell script!)

