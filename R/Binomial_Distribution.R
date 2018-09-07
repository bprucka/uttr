##################################################
# Implementation of Binomial Class               #
# Kaitlin Cornwell                               #
# August 3, 2018                                 #
##################################################

###### Binomial Functions

# Creates the basic modeling information
# arguments: predictors
# returns: tibble with model equation and type

#' Make Binomial Distribution
#' 
#' `make_binom()` creates an object that allows for analysis assuming
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
#'         that it is a binomial model and will contain the right hand side of 
#'         a model formula by creating a linear combination of predictors
#'         
#' @examples
#' make_binom(1)
#' make_binom(x)
#' make_binom(x, y)
#' 
#' @export
make_binom <- function(...) {
  # the binomial distribution is fit using base::glm
  # it takes in a formula and would set family = binomial
  
  # make the list of predictors into expressions
  predictors = enexprs(...)
  
  # output the tibble
  tibble(model_type = "Binomial", model_eq = write_equation(predictors))
}

### Binomial Base class/default

# Binomial class constructor
# arguments: model_object - the object retured from make_distribution()
#            outcome - the expression of the name of the variable to be used as the outcome of the model
#            group - an expression of the name of the variable to be used to group the data
#            opt - a list of options, the object returned from model_options()
# returns: A classed object that contains the type of model to be fit, the model equation, grouping 
#          information and model_options

#' Constructor for Binomial distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Binomial. 
Binomial <- function(model_object, outcome, group = NA, opt = NA) {
  
  # turn outcome and group variable names into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)
  
  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq))) 
  
  if (!is.na(opt))
    values <- bind_cols(values, as_tibble(opt))
  
  attr(values, "class") <- "Binomial"
} 

# Binomial fit_object()
# should not run
fit_object.Binomial <- function(object, data) {
  stop("You have not defined the type of model fit properly")
}


# Binomial model_prediction()
# should not run
#' @export
model_prediction.Binomial <- function(object, new_values = NULL) {
  stop("You have not defined the type of model fit properly")
}

# Binomial inv_transformation()
# should not run
#' @export 
inv_transformation.Binomial <- function(object, model_results, id, values) {
  stop("You have not defined the type of model fit properly")
}

#' @export
simulate_distribution.Binomial <- function(object, model_results, values, nsim, seed) {
  stop("You have not defined the type of model fit properly")
}

# Binomial make_likelihood()
# should not run
#' @export
make_likelihood.Binomial <- function(object) {
  stop("You have not defined the type of model fit properly")
}

#' @export
model_distribution.Binomial <- function(model_object, model_results, id, hist = FALSE) {
  stop("You have not defined the type of model fit properly")
}

### Binomial Frequentist class

#' Constructor for Binomial Frequentist distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Binomial.Frequestist.
#' 
#' @export 
Binomial.Frequentist <- function(model_object, outcome, group = NA, opt = NA) {
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
  structure(values, class = c("Binomial.Frequentist", "Binomial"))
}

# Fits the Binomial model with a Frequentist framework
# arguments: model_object - an object of class Binomial.Frequentist
#            model_data - the dataset to be used to fit the model
# returns: the S3 object from glm along with relevant model information including type of fit and equation
#' @export
fit_object.Binomial.Frequentist <- function(model_object, model_data) {
  
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
              do(model_results = as_result(glm(formula = !!mod_formula, family = binomial, data = .$data),
                                      model_data = .$data)) %>%
              mutate(model_type = "Binomial", fit_type = "Frequentist", model_eq = model_object$equation)
  
  if (!is.na(model_object$group))
    results <- bind_cols(groups, results)
  
  #results <- model_data %>%
  #            do(model_results = glm(formula = !!mod_formula, family = binomial, data = .)) %>%
  #            mutate(model_type = "Binomial", fit_type = "Frequentist", model_eq = model_object$equation)
  results  
}

# Predicts values for the model
# arguments: model_object - an object of type Binomial.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            id - an id value to be appended to final dataset
#            new_values - the values at which to carry out prediction at
# returns: a tibble with each row a predicted value of the dataset
#' @export
model_prediction.Binomial.Frequentist <- function(model_object, model_results, id, new_values = NULL) {
  
  if (is.na(new_values) %>% all()) {
    # get the predictors of the model
    preds <- get_predictors(model_object$equation)[-c(1,2)]
    
    # if intercept only model
    if ((preds == "1")  %>% all())
      new_values <- tibble(intercept = rep(model_results$result$coefficients, times = nrow(model_results$result$data)))
    else
      # subset data to be appended to prediction results
      new_values <- model_results$result$data %>% select(!!preds)
  }
  
  # using predict.glm function
  predict(model_results$result, new_values) %>% 
    as_tibble() %>% 
    bind_cols(new_values) %>%
    # add id column
    mutate(id = id) %>%
    rename(!!model_object$outcome[[1]] := value)
  
}

# Performs the inverse transformation of the link function
# arguments: model_object - an object of type Binomial.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            id - an id value to be appended to final dataset
#            values - the values at which to carry out the transformation of
# returns: a tibble with each row the average value at the given predictors from the model
#' @export
inv_transformation.Binomial.Frequentist <- function(model_object, model_results, id, values) {
  
  # if using success/failure parameterization
  if (grepl("cbind", colnames(model_results$model)) %>% any())
    # get the model equation and find n from it
    n <- apply(model_results$model %>% 
                  select(!!colnames(.)) %>%
                  select(starts_with("cbind")),
               1,
               sum)[1]
  # if using 0/1 parameterization
  else
    n <- 1
  
  # perform the inverse transformation
  values %>% mutate(!!model_object$outcome[[1]] := inv.logit(!!model_object$outcome[[1]]) * n)  
}

# Simulate data from a given model
# arguments: model_objects - an object of type Binomial.Frequentist
#            model_results - an S3 object of class glm contianing the fit model information
#            values - the values at which to carry out the simulation at
#            nsim - number of datasets to simulate
#            seed - value set for random number generator
# returns: a tibble with each row a simulated value at the given predictors
#' @export
simulate_distribution.Binomial.Frequentist <- function(model_objects, model_results, values, nsim = 1, seed = NULL) {
  
  if (grepl("cbind", colnames(model_results$model)) %>% any())
    n <- apply(model_results$model %>% 
                 select(!!colnames(.)) %>%
                 select(starts_with("cbind")),
               1,
               sum)[1]
  else
    n <- 1
  
  values <- values %>% 
    mutate(prob = inv.logit(!!model_objects$outcome[[1]]))
  
  set.seed(seed)
  values$prob %>%
    map(~ rbinom(nsim, n, .x)) %>%
    unlist() %>%
    as_tibble() %>%
    bind_cols(., 
              values %>% slice(rep(1:n(), each = nsim))) %>%
    select(-!!model_objects$outcome[[1]], -prob) %>%
    rename(!!model_objects$outcome[[1]] := value)
  
}

#' @export
model_distribution.Binomial.Frequentist <- function(model_object, model_results, id, hist = FALSE) {
  
  if (grepl("cbind", colnames(model_results$result$model)) %>% any())
    n <- apply(model_results$result$model %>% 
                 select(!!colnames(.)) %>%
                 select(starts_with("cbind")),
               1,
               sum)[1]
  else
    n <- 1
  
  bounds <- model_results$data %>% select(!!model_object$outcome[[1]]) %>% summarise(min = min(.), max = max(.))
  
  graphing_tbl <- tibble(x = bounds$min:bounds$max,
                         y = dbinom(x, size = n, prob = inv.logit(coef(model_results$result)[[1]])),
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
    ggtitle(glue("{id} - Binomial")) +
    theme(legend.position = "none") +
    labs(y = "Density")
}

### Binomial Bayesian Class

#' Constructor for Binomial Bayesian distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Binomial.Bayesian. 
Binomial.Bayesian <- function(model_object, outcome, group = NA, opt = NA) {
  
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)

  # get the model options specific for Bayesian analysis
  ifelse("adapt" %in% names(opt), 
         adapt <- opt$adapt,
         adapt <- 1000)
  ifelse("burn" %in% names(opt),
         burn <- opt$burn,
         burn <- 100)
  ifelse("iter" %in% names(opt),
         iter <- opt$iter,
         iter <- 1000)
  ifelse("chains" %in% names(opt),
         chains <- opt$chains,
         chains <- 3)
  ifelse("thin" %in% names(opt),
         thin <- opt$thin,
         thin <- 1)
  
  # create object
  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq)),
                   adapt = adapt,
                   burn = burn,
                   iter = iter,
                   chains = chains,
                   thin = thin) 
  if (!is.null(model_object$priors))
    values <- values %>%
                bind_cols(., model_object$priors %>% nest()) %>%
                rename(priors = data)
  
  # set the class
  structure(values, class = c("Binomial.Bayesian", "Binomial"))
}

# Fits the Binomial model with a Bayesian framework
# arguments: model_object - an object of class Binomial.Bayesian
#            model_data - the dataset to be used to fit the model
# returns: the S3 object from glm along with relevant model information including type of fit and equation
#' @export
fit_object.Binomial.Bayesian <- function(model_object, model_data) {
  
  # write the JAGS code needed to fit the model
  model_text <- create_jags_code(model_object)
  
  # Subset the data based on the model
  model_preds <- get_predictors(model_object$equation[[1]])
  
  # if there is a grouping variable given then group the data
  if (!is.na(model_object$group)) {
    # get groups to append to results dataset 
    grouping <- model_data %>% select(!!model_object$group[[1]]) %>% distinct()
    
    model_data <- model_data %>% 
                  group_by(!!model_object$group[[1]]) %>% 
                  select(!!model_object$outcome[[1]], model_preds[-1], !!model_object$group[[1]]) %>%
                  drop_na %>% 
                  nest()
  }
  else {
    # subset data
    model_data <- model_data %>% select(!!model_object$outcome[[1]], model_preds[-1]) %>% drop_na %>% nest
  }
  
  # create the model and run adapt phase
  jags_results <- model_data %>% 
                      rowwise() %>% 
                      do(results = jags.model(textConnection(model_text),
                      # the data must be given as a list
                      data = append(list(n = nrow(.$data)), as.list(.$data)),
                      n.chains = model_object$chains[[1]],
                      n.adapt = model_object$adapt[[1]],
                      inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 5)))
  
  # run the burn-in phase
  jags_results %>% rowwise() %>% do(update = update(.$results, n.iter = model_object$burn[[1]]))
  
  # get the samples and append relevant model information
  jags_results <- jags_results %>% 
                  bind_cols(model_data) %>%
                  rowwise() %>% 
                  do(model_results = as_result(coda.samples(.$results, 
                     variable.names = c(glue("beta_{model_preds}")),
                     n.iter = model_object$iter, 
                     thin = model_object$thin[[1]]),
                     model_data = .$data)) %>%
                  mutate(model_type = "Binomial", fit_type = "Bayesian", model_eq = model_object$equation)
  
  # if there was grouping
  if (!is.na(model_object$group))
    # add the column specifying groups
    jags_results <- jags_results %>% bind_cols(grouping, .) 

  jags_results
}

# Make the JAGS code for the likelihood
# arguments: an object of class Binomial.Bayesian
# returns: a string to be used as the likelihood specification in the JAGS model text
#' @export
make_likelihood.Binomial.Bayesian <- function(model_object) {
  
  # Make strings for the outcome
  outcome <- quo_name(model_object$outcome[[1]])
  outcome_string <- glue("{outcome}[i] ~ dbern(q[i])")
  
  # Make strings for the model equation
  predictors <- get_predictors(model_object$equation[[1]])
  # Start with an intercept
  predictors_string <- c("beta_int")
  # Add the other predictors plus thier values
  if (length(predictors) > 1) {
    for (i in 2:length(predictors)) {
      predictors_string <- c(predictors_string, 
                             glue("beta_{predictors[i]}*{predictors[i]}[i]"))
    }
  }
  
  predictors_string <- paste(predictors_string, collapse = " + ")
  
  logit_string <- glue("logit(q[i]) <- {predictors_string}")
    
  # Put the two pieces together
  glue("{outcome_string} \n {logit_string}")
}

#' @export
model_prediction.Binomial.Bayesian <- function(model_object, model_results, id, new_values = NA) {
  
  if (is.na(new_values) %>% all()) {
    # if intercept only model, let data = 0
    if (ncol(model_results$data) == 1) {
      pred_values <- tibble(int = rep(1, nrow(model_results$data)))
      new_values <- pred_values
    }
    else {
      pred_values <- model_results$data %>% 
                    mutate(int = 1) %>%
                    select(get_predictors(model_object$equation)[-2])
    new_values <- pred_values %>% select(-int)
    }
  }
  else {
    if (get_predictors(model_object$equation)[3] == "1")
      pred_values <- tibble(int = rep(1, nrow(new_values)))
    else
     pred_values <- new_values %>% 
        mutate(int = 1) %>%
        select(get_predictors(model_object$equation)[-2])
  }
  
  results_tbl <- model_results$result %>% 
    lapply(as.data.frame) %>% 
    reduce(bind_rows) %>% 
    as_tibble() %>%
    `colnames<-`(sub("beta_", "", colnames(.)))
  
  if (ncol(results_tbl) == 1) {
    predictions <- pred_values %>% 
                   rowwise() %>% 
                   do(!!model_object$outcome[[1]] := apply(results_tbl, 1, `*`, unlist(.)) %>%
                                                     mean)
  }
  else {
    results_tbl <- results_tbl %>%
      select(get_predictors(model_object$equation)[-2])
    
    predictions <- pred_values %>%
                   rowwise() %>%
                   do(!!model_object$outcome[[1]] := apply(results_tbl, 1, `*`, unlist(.)) %>%
                                                     apply(. , 1, mean) %>%
                                                     sum)
  }

  predictions %>%
    unnest %>%
    bind_cols(new_values)
}

#' @export
inv_transformation.Binomial.Bayesian <- function(model_object, model_results, id, values) {
  
  values %>% mutate(!!model_object$outcome[[1]] := inv.logit(!!model_object$outcome[[1]]), id = id)
  
}

### Binomial Machine Learning class

#' Constructor for Binomial Randomforest distribution
#' 
#' @param model_object an object coming from the `fit_model()` function
#' @param outcome the variable to be used as the dependent variable provided
#'        as an expression
#' @param group the variable to be used as the grouping variable provided
#'        as an expression
#' @param opt the model options provided as a function call to `model_options()`
#' 
#' @return An object of class Binomial.Randomforest
Binomial.Randomforest <- function(model_object, outcome, group = NA, opt = NA) {
  # make arguments into expressions
  outcome = enexpr(outcome)
  group = enexpr(group)
  
  # set options
  # regression or classification
  ifelse("type" %in% names(opt),
         type <- expr(!!opt$type),
         type <- expr(regression))
  # random number generator seed
  ifelse("seed" %in% names(opt),
         seed <- opt$seed,
         seed <- NA)
  # number of trees to grow
  ifelse("ntree" %in% names(opt),
         ntree <- opt$ntree,
         ntree <- 500)
  # number of variables to randomly select at each node
  ifelse("mtry" %in% names(opt),
         mtry <- opt$mtry,
         mtry <- ifelse(type == "regression", 
                        max(floor(length(get_predictors(model_object$model_eq))/3)),
                        floor(sqrt(length(get_predictors(model_object$model_eq))))))
  # replacement
  ifelse("replace" %in% names(opt),
         replace <- opt$replace,
         replace <- TRUE)
  # maximum number of terminal nodes
  ifelse("maxnodes" %in% names(opt),
         maxnodes <- opt$maxnodes,
         maxnodes <- NA)
  # should importance be tracked
  ifelse("importance" %in% names(opt),
         importance <- opt$importance,
         importance <- FALSE)

  
  
  values <- tibble(outcome = c(expr(!!outcome)),
                   group = c(expr(!!group)),
                   equation = c(expr(!!model_object$model_eq)),
                   type = c(expr(!!type)),
                   seed = seed,
                   ntree = ntree,
                   mtry = mtry,
                   replace = replace,
                   maxnodes = maxnodes,
                   importance = importance) 
  
  # set the class
  structure(values, class = c("Binomial.Randomforest", "Binomial"))
}

# Fits the Binomial model with a random forest framework
# arguments: model_object - an object of class Binomial.Randomforest
#            model_data - the dataset to be used to fit the model
# returns: the S3 object from randomForest along with relevant model information including type of fit and equation
#' @export
fit_object.Binomial.Randomforest <- function(model_object, model_data) {
  
  # check if intercept only model
  if (model_object$equation[[1]] == '1') {
    stop("You cannot fit an intercept only model using random forest")
  }
  
  # Subset the data based on the model
  model_preds <- get_predictors(model_object$equation[[1]])
  outcome <- expr(!!model_object$outcome[[1]])
  
  # if fitting classification model (need to add as.factor to outcome in formula)
  if (("type" %in% names(model_object)) && (model_object$type == "classification")) {
    # if fitting an intercept only model
    if (quo_name(model_object$equation[[1]]) == "1") 
      mod_formula <- expr(as.factor(!!outcome) ~ 1)
    else
      # if not intercept only use the given predictor
      mod_formula <- expr(as.factor(!!outcome) ~ !!model_object$equation[[1]])
  }
  # if fitting a regression model 
  else {
    if (quo_name(model_object$equation[[1]]) == "1")
        mod_formula <- expr(as.factor(!!outcome) ~ 1)
    else
      mod_formula <- expr(!!outcome ~ !!model_object$equation[[1]])
  }
  
  # if there is a grouping variable given then group the data
  if (!is.na(model_object$group)) {
    # get groups to append to results dataset
    grouping <- model_data %>% select(!!model_object$group[[1]]) %>% distinct()
    
    model_data <- model_data %>% 
      group_by(!!model_object$group[[1]]) %>% 
      select(!!model_object$outcome[[1]], model_preds[-1], !!model_object$group[[1]]) %>%
      drop_na %>%
      nest
  }
  else {
    # subset data
    model_data <- model_data %>% select(!!model_object$outcome[[1]], model_preds[-1]) %>% drop_na %>% nest
  }
  
  # get the randomforest output per group and add variables for identification
  if (!is.na(model_object$seed))
    set.seed(model_object$seed)
  # if the number of maxnodes is missing, set it equal to null
  if (is.na(model_object$maxnodes))
    maxnodes <- NULL
  else
    maxnodes <- model_object$maxnodes
  
  # fit the random forest
  results <- model_data %>%
    rowwise() %>%
    do(model_results = as_result(randomForest(formula = !!mod_formula, data = .$data, 
                                    ntree = model_object$ntree, mtry = model_object$mtry, replace = model_object$replace,
                                    maxnodes = maxnodes, importance = model_object$importance),
                                 model_data = .$data)) %>%
    mutate(model_type = "Binomial", fit_type = "Randomforest", model_eq = model_object$equation)
  
  if (!is.na(model_object$group))
    results <- results %>% bind_cols(grouping, .)
  
  results  
}

# Predicts values for the model
# arguments: model_object - an object of type Binomial.Randomforest
#            model_results - an S3 object of class glm contianing the fit model information
#            id - an id value to be appended to final dataset
#            new_values - the values at which to carry out prediction at
# returns: a tibble with each row a predicted value of the dataset
#' @export
model_prediction.Binomial.Randomforest <- function(model_object, model_results, id, new_values = NA) {

  # using predict.randomForest
  if (is.na(new_values) %>% all()) {
    predict(model_results$result) %>% 
        as_tibble() %>% 
        bind_cols(model_results$data %>% select(-heart_disease)) %>%
        # add id column
        mutate(id = id, value = as.numeric(as.character(value))) %>% 
        rename(!!model_object$outcome[[1]] := value)
  }
  else {
    predict(model_results$result, new_values) %>% 
      as_tibble() %>% 
      bind_cols(new_values) %>%
      # add id column
      mutate(id = id, value = as.numeric(as.character(value))) %>%
      rename(!!model_object$outcome[[1]] := value)
  }
  
}

# Performs the inverse transformation of the link function (already done for machine learning classification)
# arguments: model_object - an object of type Binomial.Randomforest
#            model_results - an S3 object of class randomforest contianing the fit model information
#            id - an id value to be appended to final dataset
#            values - the values at which to carry out the transformation of
# returns: a tibble with each row the average value at the given predictors from the model
#' @export
inv_transformation.Binomial.Randomforest <- function(model_object, model_results, id, values) {
  
  # return the data set since no transformation is necessary
  values
}
