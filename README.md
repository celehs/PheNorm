PheNorm
================

## Overview

The PheNorm R package provides an unsupervised phenotyping algorithm,
for electronic health record (EHR) data. A human-annotated training set
with gold-standard disease status labels is usually required to build an
algorithm for phenotyping based on a set of predictive features.
PheNorm, however, does not require expert-labeled samples for training.

The algorithm combines the most predictive variables, such as the counts
of the main International Classification of Diseases (ICD) codes, with
other EHR features. Those include for example health utilization and
processed clinical note data. PheNorm aims to obtain a score for
accurate risk prediction and disease classification. In particular, it
normalizes the surrogate to resemble gaussian mixture and leverages the
remaining features through random corruption denoising. PheNorm
automatically generates phenotyping algorithms and demonstrates the
capacity for EHR-driven annotations to scale to the next level
phenotypic big data.

The data consists of ICD codes and additional features.

The output is:

  - the predicted probability of the risk of having the phenotype

  - the coefficient beta corresponding to all the features additional to
    the ICD codes.

The main steps of the algorithm are presented in the following
flowchart:

![](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6251688/bin/ocx111f1.jpg)

### Installation

The PheNorm package can be installed using the remotes package. The
following code executed in R will get you started:

``` r
install.packages("remotes",repos = "http://cran.us.r-project.org")
```

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/g6/tcrpz9115rbfkpj7vtp_sx5w0000gn/T//RtmpjzmtNY/downloaded_packages

``` r
remotes::install_github("celehs/PheNorm")
```

    ## Skipping install of 'PheNorm' from a github remote, the SHA1 (ce356924) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library(PheNorm)
```

### Example on simulated dataset

Next, we propose a simple example in which we fit PheNorm to a simulated
dataset.

``` r
set.seed(1234)
fit.dat <- read.csv("https://raw.githubusercontent.com/celehs/PheNorm/master/data-raw/data.csv")
```

Apply the PheNorm
algorithm

``` r
fit.phenorm=PheNorm.Prob("ICD", "utl", fit.dat, nm.X = NULL, corrupt.rate=0.3, train.size=nrow(fit.dat))
head(fit.phenorm$probs)
```

    ## [1] 0.4662471 0.5384967 0.5455023 0.5419286 0.6086979 0.5320170
