identify_model <- function(data){
  cat("Step 1: Identify the statistical model generating the data\n")
  cat("Possible models: Normal, Poisson, Binomial, Exponential, etc.\n")
  cat("Output: Model object for further inference steps\n")
}

exp_check <- function(model){
  cat("Step 2: Check whether the model belongs to the Exponential Family\n")
  cat("Tasks performed:\n")
  cat("- Express pdf/pmf in canonical exponential family form\n")
  cat("- Identify natural parameter\n")
  cat("- Check full rank condition\n")
}
