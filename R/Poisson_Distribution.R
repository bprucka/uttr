##################################################
# Implementation of Poisson Class                #
# Kaitlin Cornwell                               #
# August 3, 2018                                 #
##################################################

###### Poisson

#' Make Poisson Distribution
#' 
#' `make_poisson()` creates an object that allows for analysis assuming
#' a binomial distribution.
#' 
#' The returned tibble can be row bound with other objects from the
#' `make_distribution()` functions. Then multiple models can be worked
#' with at once.
#' 
#' This tibble is the input for `set_priors()` and `fit_model()` functions.
#' 
#' @param ... A comma separated list of predictors that will be used in the
#'        model, given as unquoted expressions. 1 if fitting an intercept
#'        only model
#'
#' @return A tibble with one row. The tibble will contain information stating
#'         that it is a poisson model and will contain the right hand side of 
#'         a model formula by creating a linear combination of predictors
#'         
#' @examples
#' make_pois(1)
#' make_pois(x)
#' make_pois(x, y)
#' 
#' @export
make_pois <- function(...) {
  # the poisson distribution is fit using base::glm
  # it takes in a formula and would set family = poisson
  
  # make the list of predictors into expressions
  predictors = enexprs(...)
  
  # output the tibble
  tibble(model_type = "Poisson", model_eq = write_equation(predictors))
}

### Poisson base/default class

#' Constructor for Poisson distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Poisson 
Poisson <- function(model_object, outcome, group = NA, opt = NAA) {
  
  # turn outcome and group variable names into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)
  
  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq))) 
  
  # if there are model options then put them in the object
  if (!is.na(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  # set the class
  attr(values, "class") <- "Poisson"
}

# Poission fit_object()
# should not run
#' @export
fit_object.Poisson <- function(model_object, model_data) {
  stop("You have not defined the type of model fit properly")
}

# Poisson model_prediction()
# should not run
#' @export
model_prediction.Poisson <- function(object, new_values = NULL) {
  stop("You have not defined the type of model fit properly")
}

#' @export
simulate_distribution.Poisson <- function(model_object, model_results, values, nsim = 1, seed = NULL) {
  stop("You have not defined the type of model fit properly")
}

# Poisson inv_transformation()
# should not run
#' @export
inv_transformation.Poisson <- function(object, model_results, id, predictions) {
  stop("You have not defined the type of model fit properly")
}

#' @export
model_distribution.Poisson <- function(model_object, model_results, id, hist = FALSE) {
  stop("You have not defined the type of model fit properly")
}

### Poisson frequentist class

#' Constructor for Poisson Frequentist distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Poisson.Frequentist 
Poisson.Frequentist <- function(model_object, outcome, group = NA, opt = NA) {
  
  # turn the outcome and group variable names into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq)))
  
  # if there are model options then append them to the object
  if (!is.na(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  # set the class
  structure(values, class = c("Poisson.Frequentist", "Poisson"))
}


# Fits the Poisson model with a Frequentist framework
# arguments: model_object - an object of class Poisson.Frequntist
#            model_data - the dataset to be used to fit the model
# returns: the S3 object from glm along with relevant model information including type of fit and equation
#' @export
fit_object.Poisson.Frequentist <- function(model_object, model_data) {

  # if it is intercept only
  if (quo_name(model_object$equation[[1]]) == "1")
    # create intercept only equation
    mod_formula <- expr(!!model_object$outcome[[1]] ~ 1)
  # if it is not intercept only
  else
    # create the formula putting together the outcome and linear combination of predictors
    mod_formula <- expr(!!model_object$outcome[[1]] ~ !!model_object$equation[[1]])
  
  # if there is a grouping variable then group the data
  if (!is.na(model_object$group)) {
    groups <- model_data %>% select(!!model_object$group[[1]]) %>% distinct()
    model_data <- model_data %>% group_by(!!model_object$group[[1]]) %>% nest()
  }
  else {
    model_data <- model_data %>% nest()
  }
  
  # fit the model using glm
  results <- model_data %>%
    rowwise() %>%
    do(model_results = as_result(glm(formula = !!mod_formula, family = poisson, data = .$data),
                                 model_data = .$data)) %>%
    mutate(model_type = "Poisson", fit_type = "Frequentist", model_eq = model_object$equation)
  
  if (!is.na(model_object$group))
    results <- bind_cols(groups, results)
  
  results  
}

# Predicts values for the model
# arguments: model_object - an object of type Poisson.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            id - an id value to be appended to final dataset
#            pred_values - the values at which to carry out prediction at
# returns: a tibble with each row a predicted value of the dataset
#' @export
model_prediction.Poisson.Frequentist <- function(model_object, model_results, id, new_values = NULL) {
  
  if (is.na(new_values) %>% all()) {
    # get the predictors of the model
    preds <- get_predictors(model_object$equation)[-c(1,2)]
    
    # if intercept only model
    if ((preds == "1") %>% all()) 
      new_values <- tibble(intercept = rep(model_results$coefficients, times = nrow(model_results$data)))
    else
      # subset data to be appended to prediction results
      new_values <- model_results$data %>% select(!!preds)
  }
  
  # using predict.glm function
  predict(model_results, new_values) %>%
    as_tibble() %>%
    bind_cols(new_values) %>%
    # add id column
    mutate(id = id) %>%
    rename(!!model_object$outcome[[1]] := value)
}

# Carries out the inverse transformation of the link function
# arguments: model_object - an object of type Poisson.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            id - an id value to be appended to final dataset
#            values - the values at which to carry out the transformation of
# returns: a tibble with each row the average response at the given predictors
#' @export
inv_transformation.Poisson.Frequentist <- function(model_object, model_results, id, values) {
  
  # carry out the inverse transformation
  values %>% mutate(!!model_object$outcome[[1]] := exp(!!model_object$outcome[[1]]))  

}

# Simulate data from a fit distribution
# arguments: model_object - an object of type Poisson.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            values - the values at which to carry out simulation at
#            nsim - number of datasets to simulate
#            seed - value set for random number generator
# returns: a tibble with each row a simulated value of the dataset
#' @export
simulate_distribution.Poisson.Frequentist <- function(model_objects, model_results, values, nsim = 1, seed = NULL) {
  
  # transform the outcome into the parameter for poisson distribution
  values <- values %>% mutate(lambda = exp(!!model_objects$outcome[[1]]))
  
  # carry out the simulation starting with the predicted values
  set.seed(seed)
  values$lambda %>%
    map(~ rpois(nsim, .)) %>%
    unlist() %>%
    as_tibble() %>%
    bind_cols(., 
              values %>% slice(rep(1:n(), each = nsim))) %>%
    select(-!!model_objects$outcome[[1]], -lambda) %>%
    rename(!!model_objects$outcome[[1]] := value)
  
}

#' @export
model_distribution.Poisson.Frequentist <- function(model_object, model_results, id, hist = FALSE) {
  
  bounds <- model_results$data %>% select(!!model_object$outcome[[1]]) %>% summarise(min = min(.), max = max(.))
  
  graphing_tbl <- tibble(x = bounds$min:bounds$max,
                         y = dpois(x, lambda = exp(coef(model_results$result)[[1]])),
                         id = id)
  
  if (hist) {
    hist_data <- model_results$data %>% 
      select(!!model_object$outcome[[1]]) %>%
      group_by(!!model_object$outcome[[1]]) %>%
      summarise(n = n()) %>%
      rename(x = y)
    
    graphing_tbl <- inner_join(graphing_tbl, hist_data, by = "x")
  }
  
  model_plot <- ggplot(graphing_tbl, aes(x = x))
  
  if (hist)
    model_plot <- model_plot + geom_bar(aes(y = n/sum(n)), stat = "identity")
  
  model_plot +
    geom_line(aes(y = y), col = id) +
    ggtitle(glue("{id} - Poisson")) +
    theme(legend.position = "none") +
    labs(y = "Density")
}
