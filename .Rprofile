## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded

## also change the .libPaths and R_LIBS_USER to the right thing
.libPaths("C:/Users/zool2541/R/win-library/4.0/")

if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}
