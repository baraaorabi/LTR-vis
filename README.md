# LTR-vis
Long transcriptomic reads visualizer

## Pre-reqs

[Download](https://www.anaconda.com/distribution/) and [install](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) `conda`. 
You don't any admin/root privileges! 

Now, install `Snakemake` using `conda`
```
conda install -c bioconda snakemake
```

Clone the repo:
```
git clone https://github.com/baraaorabi/LTR-vis.git
```

## Run stuff

Check which jobs will run:
```
snakemake -np
```

Run all jobs:
```
snakemake --use-conda
```
