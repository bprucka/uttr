##################################################
# Implementation of Beta Binomial Class          #
# Kaitlin Cornwell                               #
# August 3, 2018                                 #
##################################################

###### Beta Binomial Functions

#' Make Beta Binomial Distribution
#' 
#' `make_betabinom()` creates an object that allows for analysis assuming
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
#'         that it is a beta binomial model and will contain the right hand side of 
#'         a model formula by creating a linear combination of predictors
#'         
#' @examples
#' make_betabinom(1)
#' make_betabinom(x)
#' make_betabinom(x, y)
#' 
#' @export
make_betabinom <- function(...) {
  # the binomial distribution is fit using VGAM::vglm
  # it takes in a formula and would set family = betabinom
  # make the list of predictors into expressions
  predictors = enexprs(...)
  
  # output the tibble
  tibble(model_type = "BetaBinomial", model_eq = write_equation(predictors))
}

### Beta Binomial Base/default Class

#' Constructor for Beta Binomial distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class BetaBinomial 
BetaBinomial <- function(model_object, outcome, group = NA, opt = NA) {
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq)))
  
  # if there are model options put them in the object
  if (!is.an(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  # set the class
  attr(values, "class") <- "BetaBinomial"
  values
} 

# Beta Binomial fit_object
# should not run
#' @export
fit_object.BetaBinomial <- function(model_object, model_data) {
  stop("You have not defined the type of model fit properly")
}

#' @export
simulate_distribution.BetaBinomial <- function(object, model_results, values, nsim, seed) {
  stop("You have not defined the type of model fit properly")
}

# Beta Binomial model_prediction()
# should not run
#' @export
model_prediction.BetaBinomial <- function(model_object, model_results, id, new_values = NA) {
  stop("You have not defined the type of model fit properly")
}

# Beta Binomial inv_transformation()
# should not run
#' @export
inv_transformation.BetaBinomial <- function(model_object, model_results, id, predictions) {
  stop("You have not defined the type of model fit properly")
}

#' @export
model_distribution.BetaBinomial <- function(model_object, model_results, id, hist = FALSE) {
  stop("You have not defined the type of model fit properly")
}

### Beta Binomial Frequentist class

#' Constructor for BetaBinomial Frequentist distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class BetaBinomial.Frequentist 
BetaBinomial.Frequentist <- function(model_object, outcome, group = NA, opt = NA) {
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq)))
  
  # add max option to be used if given a vector of successes
  if ("max" %in% names(opt)) {
    max <- expr(!!opt$max)
    values <- bind_cols(values, max = opt$max)
  }
    
  # set the class
  attr(values, "class") <- c("BetaBinomial.Frequentist", "BetaBinomial")
  values
}

# Fits the Beta Binomial model with a Frequentist framework
# arguments: model_object - an object of class BetaBinomial.Frequentist
#            model_data - the dataset to be used to fit the model
# returns: the S4 object from vglm along with relevant model information including type of fit and equation
#' @export
fit_object.BetaBinomial.Frequentist <- function(model_object, model_data) {
  
  # if fitting a model with success/failure outcome instead of 0/1
  if ("max" %in% names(model_object))
    outcome <- expr(cbind(!!model_object$outcome[[1]], !!model_object$max - !!model_object$outcome[[1]]))
  else
    outcome <- expr(!!model_object$outcome[[1]])
  
  # if fitting an intercept only model
  if (quo_name(model_object$equation[[1]]) == "1") 
    mod_formula <- expr(!!outcome ~ 1)
  else
    # if not intercept only use the given predictor
    mod_formula <- expr(!!outcome ~ !!model_object$equation[[1]])
  
  # if there is a grouping variable given then group the data
  if (!is.na(model_object$group)) {
    groups <- model_data %>% select(!!model_object$group[[1]]) %>% distinct()
    model_data <- model_data %>% group_by(!!model_object$group[[1]]) %>% nest()
  }
  else {
    model_data <- model_data %>% nest()
  }
  
  # get the glm per group and add variables for identification
  results <- model_data %>%
    rowwise() %>%
    do(model_results = as_result(vglm(formula = !!mod_formula, family = betabinomial, data = .$data),
                                 model_data = .$data)) %>%
    mutate(model_type = "BetaBinomial", fit_type = "Frequentist", model_eq = model_object$equation)
  
  if (!is.na(model_object$group))
    results <- bind_cols(groups, results)
  
  results  
}

# Predicts values for the model
# arguments: model_object - an object of type BetaBinomial.Frequentist
#            model_results - an S4 object of class vglm contianing the fit model information
#            id - an id value to be appended to final dataset
#            pred_values - the values at which to carry out prediction at
# returns: a tibble with each row a predicted value of the dataset
#' @export
model_prediction.BetaBinomial.Frequentist <- function(model_object, model_results, id, new_values = NA) {
  
  if (is.na(new_values) %>% all()) {
    # get the predictors of the model
    preds <- get_predictors(model_object$equation)[-c(1,2)]
    
    # if intercept only model
    if ((preds == "1") %>% all()) 
      new_values <- tibble(intercept = rep(model_results$result@coefficients[1], times = nrow(model_results$result@x)))
    else
      # subset data to be appended to prediction results
      new_values <- model_results$result@x %>% as_tibble() %>% select(-`(Intercept)`)
    
    # get predictions for original dataset
    predictions <- predictvglm(model_results$result)
  }
  else {
    # get predictions for new dataset
    predictions <- predictvglm(model_results$result, new_values)
  }
  
  # format the results
  predictions %>%
    as_tibble() %>%
    select(`logit(mu)`) %>%
    bind_cols(new_values) %>%
    # add id column
    mutate(id = id) %>%
    rename(!!model_object$outcome[[1]] := `logit(mu)`)

}

# Perform the inverse transformation of the link function
# arguments: model_object - an object of type BetaBinomial.Frequentist
#            model_results - an S4 object of class vglm contianing the fit model information
#            id - an id value to be appended to final dataset
#            predictions - the values at which to carry out prediction at
# returns: a tibble with each row the value of the mean at the given predictors
#' @export
inv_transformation.BetaBinomial.Frequentist <- function(model_object, model_results, id, predictions) {
  # get max number of trails for use to get mean
  # if the model was fit with success/failure parameterization 
  if (grepl("cbind", model_results@call$formula) %>% any()) {
    n <- str_extract(deparse(model_results@call$formula), "(?<=, )[:digit:]+(?= )")
    n <- as.numeric(n)
  }
  # if the model was fit with 0/1 parameterization
  else
    n <- 1
  
  # carry out the inverse transformation
  predictions %>% mutate(!!model_object$outcome[[1]] := inv.logit(!!model_object$outcome[[1]]) * n)
}

# Simuate values from the given model
# arguments: model_objects - an object of type BetaBinomial.Frequentist
#            model_results - an S4 object of class vglm contianing the fit model information
#            values - the values at which to carry out simulation at
#            nsim - number of datasets to simulate
#            seed - value to set for random number generation
# returns: a tibble with each row a simulated value at the given predictors
#' @export
simulate_distribution.BetaBinomial.Frequentist <- function(model_objects, model_results, values, nsim = 1, seed = NULL) {
  
  # if the model was fit with # of successes/failures
  if (grepl("cbind", model_results@call$formula) %>% any()) {
    # get n from the model statement
    n <- str_extract(deparse(model_results@call$formula), "(?<=, )[:digit:]+(?= )")
    n <- as.numeric(n)
  }
  # if the model was fit with 0/1
  else
    n <- 1
  
  # get the probability of success
  values <- values %>% 
    mutate(prob = inv.logit(!!model_objects$outcome[[1]]))
  
  # carry out the simulation
  set.seed(seed)
  values$prob %>%
    map(~ rbetabinom(nsim, n, .x)) %>%
    unlist() %>%
    as_tibble() %>%
    bind_cols(., 
              values %>% slice(rep(1:n(), each = nsim))) %>%
    select(-!!model_objects$outcome[[1]], -prob) %>%
    rename(!!model_objects$outcome[[1]] := value)
  
}

#' @export
model_distribution.BetaBinomial.Frequentist <- function(model_object, model_results, id, hist = FALSE) {
  # if the model was fit with # of successes/failures
  if (grepl("cbind", model_results$result@call$formula) %>% any()) {
    # get n from the model statement
    n <- str_extract(deparse(model_results$result@call$formula), "(?<=, )[:digit:]+(?= )")
    n <- as.numeric(n)
  }
  # if the model was fit with 0/1
  else
    n <- 1
  
  bounds <- model_results$result@y %>% 
    as_tibble() %>% 
    summarise(min = min(.), max = max(.))
  
  graphing_tbl <- tibble(x = 0:n,
                         y = dbetabinom(x,size=n,inv.logit(coef(model_results$result)[[1]]),
                                        inv.logit(coef(model_results$result)[[2]])),
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
    ggtitle(glue("{id} - Beta Binomial")) +
    theme(legend.position = "none") +
    labs(y = "Density")
}
