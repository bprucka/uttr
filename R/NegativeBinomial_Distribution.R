##################################################
# Implementation of Negative Binomial Class      #
# Kaitlin Cornwell                               #
# August 3, 2018                                 #
##################################################

###### Negative Binomial Functions

#' Make Negative Binomial Distribution
#' 
#' `make_negbinom()` creates an object that allows for analysis assuming
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
#'         that it is a negative binomial model and will contain the right hand side of 
#'         a model formula by creating a linear combination of predictors
#'         
#' @examples
#' make_negbinom(1)
#' make_negbinom(x)
#' make_negbinom(x, y)
#' 
#' @export
make_negbinom <- function(...) {
  # the negative binomial distribution is fit using VGAM::vglm
  # it takes in a formula and would set family = negbinomial
  
  # make the list of predictors into expressions
  predictors = enexprs(...)
  
  # output the tibble
  tibble(model_type = "NegativeBinomial", model_eq = write_equation(predictors))
}

### Negative Binomial Base class/default

#' Constructor for Negative Binomial distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class NegativeBinomial 
NegativeBinomial <- function(model_object, outcome, group = NA, opt = NA) {
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq))) 
  
  if (!is.na(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  attr(values, "class") <- "NegativeBinomial"
  values
} 

# Negative Binomial fit_object
# should not run
#' @export
fit_object.NegativeBinomial <- function(model_object, model_data) {
  stop("You have not defined the type of model fit properly")
}

#' @export
simulate_distribution.NegativeBinomial <- function(model_object, model_results, values, nsim = 1, seed = NULL) {
  stop("You have not defined the type of model fit properly")
}

# Negative Binomial model_prediction
# should not run
#' @export
model_prediction.NegativeBinomial <- function(model_object, model_results, id, new_values = NA) {
  stop("You have not defined the type of model fit properly")
}

# Negative Binomial inv_transformation
# should not run
#' @export
inv_transformation.NegativeBinomial <- function(model_objet, model_results, id, predictions) {
  stop("You have not defined the type of model fit properly")
}

#' @export
model_distribution.NegativeBinomial <- function(model_object, model_results, id, hist = FALSE) {
  stop("You have not defined the type of model fit properly")
}

### Negative Binomial Frequentist class

#' Constructor for Negative Binomial Frequentist distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class NegativeBinomial.Frequentist
NegativeBinomial.Frequentist <- function(model_object, outcome, group = NA, opt = NA) {
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq))) 
  
  # if model options are given add them to the object
  if (!is.na(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  # set the class
  structure(values, class = c("NegativeBinomial.Frequentist", "NegativeBinomial"))
}

# Fits the Negative Binomial model with a Frequentist framework
# arguments: model_object - an object of class NegativeBinomial.Frequentist
#            model_data - the dataset to be used to fit the model
# returns: the S4 object from vglm along with relevant model information including type of fit and equation
#' @export
fit_object.NegativeBinomial.Frequentist <- function(model_object, model_data) {
  
  # if fitting an intercept only model
  if (quo_name(model_object$equation[[1]]) == "1") 
    mod_formula <- expr(!!model_object$outcome[[1]] ~ 1)
  else
    # if not intercept only use the given predictor
    mod_formula <- expr(!!model_object$outcome[[1]] ~ !!model_object$equation[[1]])
  
  # if there is a grouping variable given then group the data
  if (!is.na(model_object$group)) {
    groups <- model_data %>% select(!!model_object$group[[1]]) %>% distinct()
    model_data <- model_data %>% group_by(!!model_object$group[[1]]) %>% nest()
  }
  else {
    model_data <- model_data %>% nest()
  }
  
  # the the glm per group and add variables for identification
  results <- model_data %>%
    rowwise() %>%
    do(model_results = as_result(vglm(formula = !!mod_formula, family = negbinomial, data = .$data),
                                 model_data = .$data)) %>%
    mutate(model_type = "NegativeBinomial", fit_type = "Frequentist", model_eq = model_object$equation)
  
  if (!is.na(model_object$group))
    results <- bind_cols(groups, results)
  
  results  
}


# Predicts values for the model
# arguments: model_object - an object of type NegativeBinomial.Frequentist
#            model_results - an S4 object of class vglm contianing the fit model information
#            id - an id value to be appended to final dataset
#            pred_values - the values at which to carry out prediction at
# returns: a tibble with each row a predicted value of the dataset
#' @export
model_prediction.NegativeBinomial.Frequentist <- function(model_object, model_results, id, new_values = NA) {
  
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
    select(`loge(mu)`) %>%
    bind_cols(new_values) %>%
    # add id column
    mutate(id = id) %>%
    rename(!!model_object$outcome[[1]] := `loge(mu)`)

}

#' @export
inv_transformation.NegativeBinomial.Frequentist <- function(model_object, model_results, id, predictions) {
  
  predictions %>% mutate(!!model_object$outcome[[1]] := exp(!!model_object$outcome[[1]]))

}

#' @export
simulate_distribution.NegativeBinomial.Frequentist <- function(model_objects, model_results, values, nsim = 1, seed = NULL) {
  
  model_values <- model_results@x %>% as_tibble() %>% select(-`(Intercept)`)
  values_noout <- values %>% select(-id, -!!model_objects$outcome[[1]])
  
  if (((dim(model_values) == dim(values_noout)) %>% all) && ((model_values == values_noout) %>% all)) 
    values <- values %>% mutate(mu = exp(!!model_objects$outcome[[1]]), size = exp(predictvglm(model_results)[,2]))
  else
    values <- values %>% mutate(mu = exp(!!model_objects$outcome[[1]]), size = exp(predictvglm(model_results, values)[,2]))
  
  set.seed(seed)
  values %>% 
    rowwise() %>% 
    do(simulated = rnbinom(n = nsim, size = .$size, mu = .$mu)) %>%
    unlist() %>%
    as_tibble() %>%
    bind_cols(.,
              values %>% slice(rep(1:n(), each = nsim))) %>%
    select(-!!model_objects$outcome[[1]], -mu, -size) %>% 
    rename(!!model_objects$outcome[[1]] := value)
  
}

#' @export
model_distribution.NegativeBinomial.Frequentist <- function(model_object, model_results, id, hist = FALSE) {
  
  bounds <- model_results$result@y %>% 
            as_tibble() %>% 
            summarise(min = min(.), max = max(.))
  
  graphing_tbl <- tibble(x = bounds$min:bounds$max,
                         y = dnbinom(x,size=exp(coef(model_results$result)[[2]]),
                                     mu=exp(coef(model_results$result)[[1]])),
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
    ggtitle(glue("{id} - Negative Binomial")) +
    theme(legend.position = "none") +
    labs(y = "Density")
}
