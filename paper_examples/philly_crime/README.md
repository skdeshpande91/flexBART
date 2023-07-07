## Philly crime data analysis

The file `philly_crime_data.RData` contains the pre-processed crime data.
The script `philly_loo.R` runs a single leave-one-out fold (as reported in Section 4.2 of the paper).
You must specify the heldout tract (the object `vertex` in the script) and the specific implementation (the object `method`) in the script.
Note `vertex` must be between 1 and 384 and `method` must be one of `"networkBART"`, `"flexBART"`, `"BART-default"`, or `"BART-alt"`. 
The file `helpers.R` contains a number of auxiliary functions for generating posterior predictive draws and for computing error metrics.
The file `philly_cv.R` runs a single cross-validation fold of the experiment reported in Appendix A2.
