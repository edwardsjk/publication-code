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

`helper.R`

-   R helper functions for primary analysis

`example_analysis.R`

-   R code for all analyses
