# The Rogan-Gladen estimator for outcome misclassification

## Files 

`data/rg_dat.csv`: Main study data, including variables

-   `id`: Participant identification number
-   `ystar`: Potentially misclassified HIV rapid test result
-   `r`: Indicator of available data on `ystar`
-   `drinking`: Indicator of self-reported alcohol consumption in past year
-   `edu`: Self reported education status
-   `selfreportsw`: Indicator of report of accepting cash for sex in past 12 months
-   `preg`: Indicator of pregnancy
-   `age_g30`: Indicator of age \> 30
-   `sti`: Indicator of STI symptoms in past 4 weeks

`R/helper.R`

- R helper functions for primary analysis

`R/example_analysis.R`

- analysis code using R

`python/example_analysis.py` 

- analysis code using python

## System Details

R 4.5.1

- Dependencies: `rootSolve_1.8.2.4`, `tidyr_1.3.1`, `dplyr_1.1.4 `, and `fastDummies_1.7.5`

Python 3.9.4

- Dependencies: `numpy (1.25.2)`, `pandas (1.4.1)`,  `formulaic (0.5.2)`,  `delicatessen (3.1)`

