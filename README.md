![Build check](https://github.com/morinlab/GAMBLR.viz/actions/workflows/build_check.yaml/badge.svg)
![GitHub R package version](https://img.shields.io/github/r-package/v/morinlab/GAMBLR.viz)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/morinlab/GAMBLR.viz)
![GitHub last commit](https://img.shields.io/github/last-commit/morinlab/GAMBLR.viz)


# GAMBLR.viz

Collection of functions to make plots for Genomic Analysis of Mature B-cell
Lymphomas in R.

For detailed information about how to install, bundled resources, functionality
overview, tutorials, frequently asked questions and much more please visit the
website [morinlab.github.io/GAMBLR.viz](https://morinlab.github.io/GAMBLR.viz/)

# Contributing

## Cloning repo for the code development

The easiest way to obtain and contribute to GAMBLR.viz is to do this via cloning the repository

```
cd
git clone git@github.com:morinlab/GAMBLR.viz.git
```

In your R editor of choice, set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR.viz")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install()
```

As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.

## Function conflicts

This package relies on the use of some functions (e.g. `get_gambl_metadata()`, `get_coding_ssm()` etc) that exist in 2 different versions: [GAMBLR.data](https://github.com/morinlab/GAMBLR.data) for the users who do not have access to GSC and [GAMBLR.results](https://github.com/morinlab/GAMBLR.results) for the Morin Lab users with access to GSC. If your contribution relies on the use of such functions, please follow these 2 steps:

* *DO NOT* prepend the function use with `<package>::` (for example, `<package>::function()`), and
* *DO NOT* add the corresponding package to the `@import` section of the function

Following these steps will ensure correct usage of the proper function depending on which package is loaded in the session and will avoid functionality conflicts.
