# Defines a BatchtoolsParam object for submission on a SLURM cluster.

BPPARAM <- BatchtoolsParam(10, 
    cluster="slurm", template="parallel/slurm-aaron.tmpl", 
    logdir="parallel", log=TRUE, 
    RNGseed=10000L,
    resources=list(walltime=20000, memory=8000, ncpus=1))
