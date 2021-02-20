
# PheNorm: Unsupervised Gold-Standard Label Free Phenotyping Algorithm for EHR Data

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
normalizes the surrogate to resemble Gaussian mixture and leverages the
remaining features through random corruption denoising. PheNorm
automatically generates phenotyping algorithms and demonstrates the
capacity for EHR-driven annotations to scale to the next level
phenotypic big data.

The input data consists of ICD codes and additional features.

The PheNorm output includes:

-   the predicted probability of the risk of having the phenotype

-   the coefficient beta corresponding to all the features additional to
    the ICD codes.

The main steps of the algorithm are presented in the following
flowchart:

![](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6251688/bin/ocx111f1.jpg)

## Installation

The PheNorm package can be installed from CRAN or GitHub. The following
code executed in R will get you started:

### Stable Version

Install stable version from CRAN:

``` r
install.packages("PheNorm")
```

### Development Version

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("celehs/PheNorm")
```

## Reference

Yu S, Ma Y, Gronsbell J, Cai T, Ananthakrishnan AN, Gainer VS, Churchill
SE, Szolovits P, Murphy SN, Kohane IS, Liao KP, Cai T. Enabling
phenotypic big data with PheNorm. J Am Med Inform Assoc. 2018 Jan
1;25(1):54-60. doi: 10.1093/jamia/ocx111. PMID: 29126253; PMCID:
PMC6251688. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6251688/>
