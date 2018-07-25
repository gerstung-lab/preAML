# Prediction of AML development in healthy individuals
This repository contains code accompanying [_Prediction of AML risk in healthy individuals_ by Abelson, Collord, et al. (2018), _Nature_, **559**:400-404. http://dx.doi.org/10.1038/s41586-018-0317-6](http://dx.doi.org/10.1038/s41586-018-0317-6).

All code is contained in `preAML.R` with the corresponding output in `preAML.html`.

## Requirements
The script was run in `R` versions `3.4.1` and `3.5.1`.
```R
install.packages(c("survAUC","survivalROC","stringr","dplyr","readr","devtools","ROCR","rmarkdown"))
devtools::install_github("gerstung-lab/CoxHD/CoxHD@bc60c16")
devtools::install_github("mg14/mg14@6a63283")
```

Necessary data files are contained in the `data` directory.

## Run
The provided output has been generated in `R-3.5.1` from the toplevel directory using the command
```R
rmarkdown::render("preAML.R")`
```

The typical run time is about 2h, using 4 cores, the longest part is section `7.4.8`, which may be omited by commenting out the corresponding sections in `preAML.R`. 

