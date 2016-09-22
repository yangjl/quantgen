# farmeR

This is an R packages to generate genomic and bioinformatic pipelines and submit jobs on HPC system running slurm.

## Installation

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `imputeR` from github.

```R
#install.packages(devtools)
devtools::install_github("yangjl/quantgen")
library(quantgen)
```

List all the functions in the package and find help.

```R
ls(getNamespace("quantgen"), all.names=TRUE)
```

## License
It is a free and open source software, licensed under [GPLv3](LICENSE).
