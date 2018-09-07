##################################################
# Complete Pipeline                              #
# Kaitlin Cornwell                               #
# August 3, 2018                                 #
##################################################

# Step 1: create an array of model types
# this is a list of functions whose result would be binded together
#   each model would take in type of model (defined by what function is called)
#   and any predictors and other options that would be needed

# Step 2: fit models to data by group
# this is where you specify what data set you are using
#   result would be a tibble with #models x #groups rows 
#   and would contain the object returned by the model

# Step 3: do any calculations or plotting
#   use (user or predefined) functions to perform other actions
#   preliminary functions to be implemented are: predict() and simulate()

##############################################################################

### Default/Base Functions


#' Make an arbitrary distribution
#' 
#' `make_distribution()` returns an object that represents a distribution
#' 
#' This is a generic function which takes in a distribution name and a set
#' comma separated predictiors and returns a tibble with the information.
#' 
#' @param dist an expression containing the name of the distribution
#' @param ... a comma separated list of predictors that will be in the
#'        model, given as unquoted expressions
#'        
#' @return A tibble with one row. The tibble will contain information stating
#'         the name of the distribution and the right hand side of the model
#'         formula by forming a linear combination of the given predictors.
#'         
#' @examples 
#' make_distribution(binomial, 1)
#' make_distribution(normal, weight, height)
#' @export
make_distribution <- function(dist, ...) {
  
  # get distribution and list of predictors as expressions
  dist <- enexpr(dist)
  predictors = enexprs(...)
  
  # output the tibble
  tibble(model_type = quo_name(dist), model_eq = write_equation(predictors))
}

# Takes in a list of predictors and returns them as a linear combination
# arguments: list of predictors as expressions
# returns: predictors as an equation
write_equation <- function(predictors) {
  
  # if we get only a 1 we must turn it into a list so it will bind to other rows of type list
  if (length(predictors) == 1 && predictors[[1]] == 1)
    return(c(as.name(predictors[[1]])))
  
  # make the equation
  model_equation <- reduce(predictors, ~expr(!!.x + !!.y))
  
  # return the equation as an atomic vector - will be coerced to a list in final output tibble
  # this is necessary since a tibble column cannot be of type symbol/language
  return(c(expr(!!model_equation)))
}

# Gets a list of the predictors in a model
# To be used internally
# arguments: a model equation as a string
# returns: a list of the predictors to be included in the model

#' Returns a vector of predictors
#' 
#' `get_predictors()` returns a list with the variable names found in the model equation.
#' 
#' @param model_eq A model equation as a string
#' 
#' @return A vector with each variable found in model_eq as a string
#' 
#' @examples 
#' get_predictors("y~1")
#' get_predictors("y~height + weight")
#' 
#' @export
get_predictors <- function(model_eq) {
  # if it is an intercept only model
  if (model_eq == "1")
    model_vars <- c("int")
  # if it is not an intercept only model
  else 
    model_vars <- str_extract_all(as.character(model_eq), boundary("word")) %>% 
      unlist() %>%
      c("int", .)
  
  return(model_vars)
}

#' Assigns predictors to a model
#' 
#' `set_priors()` sets the priors for each variable in a model.
#' 
#' Each prior is written following JAGS notation and specified using 
#'     `variable = JAGS`. User must specify a prior for the intercept 
#'     using `int=`. The priors are model specific so if the user
#'     supplies priors for variables that are not found in the model
#'     equation they will not be retained.
#'     
#' @param model_array a tibble containing the model information (distribution
#'        and model equation). Often the output of a `make_distribution()` 
#'        function.
#' @param ... a list of expressions set up as "predictor = distribution" 
#'           where the distribution is written in JAGS syntax
#'           
#' @return a modified version of model_array that has the relevant predictors
#'          for each model(row)
#'          
#' @export
set_priors <- function(model_array, ...) {
  # Make the list of priors into a list of expressions
  priors_given <- enexprs(...)
  
  # Ensure an intercept prior was defined
  if (!("int" %in% names(priors_given)))
    stop("Must define a prior for the intercept")
  
  # For each row set the priors
   model_array %>%
      rowwise() %>%
      do(prior = set_individual_priors(.$model_eq, priors_given)) %>%
      bind_cols(model_array, .) %>%
      unnest(prior) %>%
      rename(priors = data)
}

# Helper function for set_priors that creates the code to set the 
#     priors within RJAGS
# arguments: the model equation and a list of prior distributions for
#            each predictor in the equation
# returns: a tibble of the priors for each predictor
set_individual_priors <- function(model_eq, prior_list) {
  
  # Get the relevant predictors
  model_vars <- get_predictors(model_eq)
  
  prior_list <- prior_list[model_vars]
  
  # For each predictor the expression from prior_list, "pred = dist" into "pred ~ dist"
  # and return them as a tibble
  lapply(prior_list, deparse) %>%
      tibble(., name = names(prior_list)) %>%
      mutate(priors = glue("{name} ~ {.}")) %>%
      select(priors) %>%
      nest()
}

# Method dispatch for fit_object
fit_object <- function(model_object, model_data) {
  UseMethod("fit_object")
}

# Default fit_object function
#' @export
fit_object.default <- function(model_object, model_data) {
  return()
}

#' Fit the specified model
#' 
#' Fits the specified model and returns a tibble containing model fit information.
#' 
#' Currently, only the binomial distribution is supported within the Bayesian
#' and random forest frameworks. All models with a specific `make_distribution()`
#' function can be fit within a Frequentist framework. Priors must be set before
#' fitting a model within a Bayesian framework.
#' 
#' If a grouping variable is specified then each model is fit separately to each 
#' group contained within the grouping variable. 
#' 
#' @param model_array a tibble of objects retreived from the `make_distribution()` function
#' @param method the framework with which the model shoudl be fit provided as an unquoted 
#'        argument. Either `frequentist`, `bayeasian` or `randomforest`.
#' @param model_data the data frame or tibble containing the variables of the model
#' @param outcome the name of the varaible to act as the dependent variable of the model 
#'        specified as an expression
#' @param group teh name of teh variable to act as the grouping variable of the data 
#'        specified as an expression.
#' @param opt additional model options provided as a fucntion call to `model_options()`
#' 
#' @return A tibble contianing a model_results object (contains the given data set and 
#'         the model fit), along with relevant information about the model including
#'         the distribution the model was fit to, what variable was set as the dependent
#'         variable, which group (if not NA) the results correpond to, and the formula
#'         for the model.
#'         
#' @examples
#' trial_data <- tibble(w = runif(20, 0, 2), x = sample(c(1, 2), 20, TRUE),
#'                      y = rbinom(20, 1, .75), z = rbinom(20, 15, .5))
#' make_pois(1) %>% fit_model(frequentist, trial_data, z)
#' make_pois(w, x) %>% fit_model(frequentist, trial_data, z)
#' make_pois(w) %>% fit_model(frequentist, trial_data, z, x)
#' make_binom(1) %>% fit_model(frequentist, trial_data, y)
#' make_binom(w) %>% fit_model(frequentist, tiral_data, z, opt = model_options(max = 15))
#' make_binom(w) %>% set_priors(int = dnorm(0, .01), w = dunif(-100, 100)) %>%
#'                   fit_model(bayesian, trial_data, y)
#'         
#' @export
fit_model <- function(model_array, method, model_data, outcome, group = NA, opt = NA) {
  # set up an object with appropriate classes
  # define functions to match with base class and child classes
  # call the functions
  
  # make method into a string with the first letter uppercase
  method <- enexpr(method)
  method <- quo_name(method)
  method <- first_up(method)
  
  # make sure priors were defined if running a bayesian analysis
  if (method == "Bayesian" && !("priors" %in% names(model_array)))
    stop("You must define priors before running a Bayesian analysis")
  
  # make outcome and group into expressions
  outcome <- enexpr(outcome)
  group <- enexpr(group)

  # model_objects is a tibble with each element of its appropriate class
  model_objects <- model_array %>%
    # get class name to be called later
    mutate(type = paste0(.$model_type, ".", method)) %>%
    # evaluate each row separately
    rowwise() %>%
    # make class objects
    do(object = eval(call(.$type, ., outcome, group, opt)))
  
  # create the tibble that will be returned with fit results
  model_results <- model_objects %>%
    rowwise() %>%
    # fit each row
    do(results = fit_object(.$object, model_data)) %>%
    # unnest so each group's results has its own row
    unnest()
  
  #suppressWarnings(model_results %>% rowwise() %>%
  #                 # combine the outcome and linear combination of predictors
  #                mutate(model_equation = glue("{outcome}~{.$model_eq}")) %>%
  #                select(-model_eq))
  
  model_results %>% rowwise() %>%
    # combine the outcome and linear combination of predictors
    mutate(model_equation = paste0(quo_name(outcome), "~", quo_name(model_eq))) %>%
    select(-model_eq)
}

# Helper function to setting the options for the models
# arguments: a comma separated list specifiying the option and the value
# returns: a named list of the option and its value as expressions

#' Set options for a model
#' 
#' `model_options()` sets model specific options for use within `model_fit()`.
#' 
#' @param ... A list of unquoted expressions setting the options in the form 
#'        of `option = value`.
#'        
#' @return Returns a named list containing the model options specified
#' 
#' @export
model_options <- function(...) {
  # Get the arguments as a list
  opt <- enexprs(...)
  # return the list
  return(opt)
}

as_result <- function(model_result, model_data) {
  results_list <- list(result = model_result, data = model_data)
  class(results_list) <- "model_results"
  #model_list <- model_list %>% map(function(x) {class(x) <- "model_results"; return(x)})
  results_list
}

#' Predict values from a fitted distribution
#' 
#' `do_prediction()` returns the predicted value for a set of predictors for a fitted model
#' 
#' Values are returned on the original scale of the outcome variable.
#' 
#' If a group was given in `fit_model()` then a prediction is given for each group separately.
#' 
#' @param model_array a tibble retreived from the `fit_model()` function
#' @param new_values a data frame or tibble contianing the values of each predictor to predict
#'        at. Defaults to the given dataset.
#'        
#' @return A tibble containing each predicted value along with which values were used for the
#'         prediction. Additional model information is returned including the type of distribution
#'         the model equation, and grouping information.
#'         
#' @examples 
#' trial_data <- tibble(w = runif(20, 0, 2), x = sample(c(1,2), 20, TRUE), y = rbinom(20, 15, .75))
#' new_data <- tibble(w = runif(5, 1, 2.5))
#' 
#' pois_model <- make_pois(w) %>% fit_model(frequentist, trial_data, y)
#' pois_model %>% do_prediction()
#' pois_model %>% do_prediction(new_data)
#' pois_model %>% do_prediction(new_data, 3)
#' 
#' @export
do_prediction <- function(model_array, new_values = NA) {
  # if there are 4 columns then there was no grouping variable
  if (ncol(model_array) == 4) {
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type, 
                            .,
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1])))) %>%
      ungroup() %>%
      # add id column to merge with prediction results
      mutate(id = row_number())
  }
  else {
    # get grouping variable
    grouping <- colnames(model_array %>% select(-model_results, -model_type, -fit_type, -model_equation))
    grouping <- as.name(grouping)
    
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type, 
                            .,
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1]),
                            group = grouping))) %>%
      ungroup() %>% 
      # add id column to merge with prediction results
      mutate(id = row_number())
  }
  

  complete_array <- bind_cols(model_objects, model_array)
  
  # do the prediction and make each simulated value into its own row
  complete_array %>%
    rowwise() %>%
    # do the prediction (using proper method)
    do(predict = model_prediction(.$object, .$model_results, .$id, new_values)) %>%
    bind_cols(complete_array) %>%
    rowwise() %>%
    # transform the outcome to be the mean (using proper method)
    do(transformation = inv_transformation(.$object, .$model_results$result, .$id, .$predict)) %>%
    # format results
    unnest() %>%
    inner_join(., complete_array, by = "id") %>%
    select(-object, -id, -model_results)
}

# Method dispatch for model_prediction
# carries out predictions for a given model
model_prediction <- function(model_object, model_results, id, new_values = NA) {
  UseMethod('model_prediction')
}

# Default function for model_prediction
#' @export
model_prediction.default <- function(model_object, model_results, id, new_values = NA) {
  return()
}

# Method dispatch for inv_transformation
# performs the inverse transformation of the link function to return the mean of the distribution
inv_transformation <- function(model_object, model_results, id, values) {
  UseMethod('inv_transformation')
}

# Default function for inv_transformation
#' @export
inv_transformation.default <- function(model_object, model_results, id, values) {
  return()
}

#' Simulate from a fitted model
#' 
#' Simulate datasets from the model previously fit. This function has the 
#' additional capability of simulating new data. 
#' 
#' This simulation function first uses the model to determine the average
#' value at a given set of predictors. Then it generates a random number
#' from the specified type of distribution using the parameters retreived
#' from the prediction.
#' 
#' If a grouping variable was indicated then the simulation is carried out
#' separately for each group.
#' 
#' This function is only supported for the frequentist methods.
#' 
#' @param model_array a tibble retreived from the `fit_model()` function
#' @param new_values a data frame or tibble containing the values of each
#'        predictor to simulate at. Defaults to the given dataset.
#' @param nsim the number of datasets to simulate
#' @param seed an integer specifying how to initialize the random number
#'        generator
#'
#' @return A tibble with each row representing one simulated value. Additional
#'         model information is returned including the type of model, the
#'         equation corresponding to the model, and teh value of the predictors
#'         at which a given value was simulated from.
#'         
#' @examples 
#' trial_data <- tibble(w = runif(20, 0, 2), x = sample(c(1,2), 20, TRUE), y = rbinom(20, 15, .75))
#' new_data <- tibble(w = runif(5, 1, 2.5))
#' 
#' pois_model <- make_pois(w) %>% fit_model(frequentist, trial_data, y)
#' pois_model %>% do_simulation()
#' pois_model %>% do_simulation(new_data)
#' pois_model %>% do_simulation(new_data, 3)
#' 
#' @export
do_simulation <- function(model_array, new_values = NA, nsim = 1, seed = NULL) {
  
  # if there are 4 columns then there was no grouping variable
  if (ncol(model_array) == 4) {
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type, 
                            .,
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1])))) %>%
      ungroup() %>%
      mutate(id = row_number())
  }
  else {
    # get grouping variable
    grouping <- colnames(model_array %>% select(-model_results, -model_type, -fit_type, -model_equation))
    grouping <- as.name(grouping)
    
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type,
                            ., 
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1]),
                            group = grouping))) %>%
      ungroup() %>%
      mutate(id = row_number())
  }
  
  # add an id column to model_array to merge with simulation results
  complete_array <- bind_cols(model_objects, model_array)
  
  complete_array %>%
      rowwise() %>%
      do(predicted = model_prediction(.$object, .$model_results, .$id, new_values)) %>%
      bind_cols(., complete_array) %>%
      do(simulated = simulate_distribution(.$object, .$model_results$result, .$predicted, nsim, seed)) %>%
      unnest() %>%
      inner_join(., complete_array, by = "id") %>%
      select(-id, -model_results, -object)
}

# Method dispatch for simulate_distribution
# simualtes new data from a distribution
simulate_distribution <- function(model_object, model_results, values, nsim, seed) {
  UseMethod('simulate_distribution')
}

# default method for simulate_distribution
#' @export
simulate_distribution.default <- function(model_object, model_results, values, nsim, seed) {
  return()
}

# Method dispatch for make_likelihood
# writes the RJAGS code for the likelihood of an object
make_likelihood <- function(model_object) {
  UseMethod('make_likelihood')
}

# Default function for make_likelihood
make_likelihood.default <- function(model_object) {
  return()
}

#' Plot each distribution
#' 
#' `plot_distribution()` plots each distribution on a separate graph and returns
#' a tibble with the graphs.
#' 
#' If a grouping variable was supplied to `fit_model()` a graph is made for 
#' each group.
#' 
#' @param model_array a tibble retreived from the `fit_model()` function
#' @param hist a logical indicating whether a histogram should overlay the
#'        distribution
#' 
#' @return A tibble with each row representing one distribution. The `graphics`
#'         column contains a ggplot object with the specific graph. Additional
#'         model information is returned including the type of model and the
#'         equation corresponding to the model.
#'         
#' @export
plot_distribution <- function(model_array, hist = FALSE) {
  # if there are 4 columns then there was no grouping variable
  if (ncol(model_array) == 4) {
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type, 
                            .,
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1])))) %>%
      ungroup() %>%
      # add id column to merge with prediction results
      mutate(id = row_number())
  }
  else {
    # get grouping variable
    grouping <- colnames(model_array %>% select(-model_results, -model_type, -fit_type, -model_equation))
    grouping <- as.name(grouping)
    
    # create classed objects
    model_objects <- model_array %>%
      ungroup() %>%
      mutate(type = paste0(.$model_type, ".", .$fit_type)) %>%
      rowwise() %>%
      do(object = eval(call(.$type, 
                            .,
                            outcome = as.name(strsplit(.$model_equation, "~")[[1]][1]),
                            group = grouping))) %>%
      ungroup() %>% 
      # add id column to merge with prediction results
      mutate(id = row_number())
  }
  
  complete_array <- bind_cols(model_objects, model_array)
  
  # create the ggplot
  complete_array %>%
    rowwise() %>%
    do(graphics = model_distribution(.$object, .$model_results, .$id, hist)) %>%
    bind_cols(model_array, .) %>%
    select(-model_results)
}

# method dispatch for model_distribution()
model_distribution <- function(model_object, model_results, id, hist = FALSE) {
  UseMethod('model_distribution')
}

# default method for model_distribution()
#' @export
model_distribution.default <- function(model_object, model_results, id, hist = FALSE) {
  return()
}

# Combine all individual distribution plots
# arguments: model_array - the object created by fit_model()
# returns: a ggplot object containing all of the distributions overlayed on one graph

#' Combine all individual distribution plots
#' 
#' `combine_plots()` combines each plot and overlay them onto one graph. 
#' 
#' If a grouping varaible was specified in `fit_model()` one graph is made for 
#' each group.
#' 
#' @param model_array a tibble retreived from the `fit_model()` function
#' @param hist A logical value indicating whether a histogram of the data
#'        should be overlayed on the graph
#' 
#' @return A ggplot object containing the overlayed graphs.
#' 
#' @export
combine_plots <- function(model_array, hist = FALSE) {
  
  # create the dataset for the ggplot
  model_array <- model_array %>%
                  rowwise() %>%
                  # get the data out of each individual ggplot object
                  do(data = use_series(., graphics) %>% use_series(., data)) %>%
                  bind_cols(model_array) %>%
                  unnest(data) 
  
  # if there is grouping change the id
  if (colnames(model_array)[1] != "model_type")
    model_array <- model_array %>% 
                  mutate(id = ifelse(id %% 2 == 0, id - 1, id)) %>%
                  mutate(id = ceiling(id/2)) 
  
  model_array <- model_array %>%
                 mutate(complete_id = glue("{.$id} - {.$model_type}")) %>%
                 mutate(complete_id= as.factor(complete_id))
  
  dist_plots <- ggplot(model_array, aes(x = x))
  
  if (hist) {
    if (!("n" %in% colnames(model_array))) {
      stop("You must first make individual graphs with histograms overlaying.
           Use hist = TRUE in plot_distribution()")
    }
    
    dist_plots <- dist_plots + geom_bar(aes(y = n/sum(n)), stat = "identity")
  }
  
  # if there is grouping need to separate graphs by group
  if (colnames(model_array)[1] != "model_type") {
    # create the ggplot
    dist_plots + 
      geom_line(aes(y = y, col = complete_id)) +
      labs(y = "Density", color = "Distribution", x = strsplit(model_array$model_equation[[1]], "~")[[1]][1]) +
      facet_wrap(colnames(model_array)[1])
  }
  else {
    # create the ggplot
    dist_plots + 
      geom_line(aes(y = y, col = complete_id)) +
      labs(y = "Density", color = "Distribution", x = strsplit(model_array$model_equation[[1]], "~")[[1]][1])
  }
  
}

### Helper Functions

# make first letter uppercase - to be used internally to create classes
# arguments: x - a string
# returns: x with the first letter as an uppercase
first_up <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Write the jags coded needed to run RJAGS
# arguments: model_object - a tibble returned by the make_distribution() and set_priors() functions
# returns: A string containing the properly written and formatted RJAGS code
create_jags_code <- function(model_object) {
  
  # Get the likelihood's RJAGS code
  like <- make_likelihood(model_object)
  
  # Get the prior's RJAGS code
  priors <- make_prior(model_object$priors[[1]])
  
  # paste the the likelihood and prior specifications together with relevant formatting
  glue("model { \n
        # Likelihood \n
        for (i in 1:n) { \n
          <<like>>
        }
        # Priors \n
        <<priors>>
      }", .open = "<<", .close = ">>")

}

# Create the RJAGS code for the prior distributions
# arguments: the priors as a tibble
# returns: a string of the priors in thier proper format for RJAGS
#          with each predictor taking the form beta_predictor
make_prior <- function(priors_tbl) {
  
  # Make it so the coefficients take on the form beta_predictor
  priors_list <- glue("beta_{unlist(priors_tbl)}")
  
  # Put all of the priors in one string separated by a new line (RJAGS format)
  paste(priors_list, collapse = " \n ")
  
}
