# Function to calculate % of variation between average fitted values at the maximum and minimum values of the predictor "x"
# It computes the average % of variation of the response caused by the predictor "x"
# Formula: ((max-min)/min)*100

# Input:
# - model: name of the model
# - min: data frame with the minimum value of "x" and the mean value for the other predictors (if any)
# - max: data frame with the maximum value of "x" and the mean value for the other predictors (if any)
# - exp: LOGICAL; to back-transform when the variable is modelled on the natural-log scale (defaults to TRUE)

# Value: the average variation

percent_variation <- function(model, min, max, exp = TRUE) {
  if (!exp) {
    ((fitted(model, max, re_formula = NA)[1]-fitted(model, min, re_formula = NA)[1])/fitted(model, min, re_formula = NA)[1])*100
  } else {
    ((exp(fitted(model, max, re_formula = NA)[1])-exp(fitted(model, min, re_formula = NA)[1]))/exp(fitted(model, min, re_formula = NA)[1]))*100
  }
}
