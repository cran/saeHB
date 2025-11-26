#' @title Small Area Estimation using Hierarchical Bayesian under Beta Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Beta Distribution. The range of data must be \eqn{0<y<1}. The data proportion is supposed to be implemented with this function.
#' @param formula Formula that describe the fitted model
#' @param ni Data column name indicating domain sizes. In out-of-sample areas, sizes must be 0.
#' @param deff_iw The factor or string deff_iw represents an approximation of the design effect for area-i, associated with the proportion of area i.
#' @param iter.update Number of updates with default \code{3}
#' @param iter.mcmc Number of total iterations per chain with default \code{10000}
#' @param coef a vector contains prior initial value of Coefficient of Regression Model for fixed effect with default vector of \code{0} with the length of the number of regression coefficients
#' @param var.coef a vector contains prior initial value of variance of Coefficient of Regression Model with default vector of \code{1} with the length of the number of regression coefficients
#' @param thin Thinning rate, must be a positive integer with default \code{2}
#' @param burn.in Number of iterations to discard at the beginning with default \code{2000}
#' @param tau.u  Prior initial value of inverse of Variance of area random effect with default \code{1}
#' @param data The data frame
#' @param seed number used to initialize a pseudorandom number generator (default seed = 1). The random number generator method used is "base::Wichmann-Hill".
#' @param quiet if TRUE, then messages generated during compilation will be suppressed (default TRUE).
#' @param plot if TRUE, the autocorrelation, trace, and density plots will be generated (default TRUE).
#'
#' @return This function returns a list of the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{plot}{Trace, Dencity, Autocorrelation Function Plot of MCMC samples}
#'
#' @export BetaDeff
#'
#' @examples
#' \donttest{
#' # Data Generation
#' set.seed(123)
#' m <- 30
#' x1 <- runif(m, 0, 1)
#' x2 <- runif(m, 0, 1)
#' x3 <- runif(m, 0, 1)
#' x4 <- runif(m, 0, 1)
#' b0 <- b1 <- b2 <- b3 <- b4 <- 0.5
#' u <- rnorm(m, 0, 1)
#' pi <- rgamma(1, 1, 0.5)
#' yi <- exp(b0 + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4 + u)
#' Mu <- yi / (1 + yi)
#' A <- Mu * pi
#' B <- (1 - Mu) * pi
#' y <- rbeta(m, A, B)
#' y <- ifelse(y < 1, y, 0.99999999)
#' y <- ifelse(y > 0, y, 0.00000001)
#' MU <- A / (A + B)
#' vardir <- A * B / ((A + B)^2 * (A + B + 1))
#' dataBeta <- as.data.frame(cbind(y, x1, x2, x3, x4, vardir))
#' dataBeta$ni <- sample(20:100, m, replace = TRUE)
#' dataBeta$deff <- round(runif(m, 1, 2), 2)
#' dataBetaNs <- dataBeta
#' dataBetaNs$y[c(3, 14, 22, 29, 30)] <- NA
#' dataBetaNs$vardir[c(3, 14, 22, 29, 30)] <- NA
#'
#'
#' ## Compute Fitted Model
#' ## y ~ x1 + x2
#'
#'
#' ## For data without any nonsampled area
#'
#' vc <- c(1, 1, 1)
#' c <- c(0, 0, 0)
#' formula <- y~x1 + x2
#' dat <- dataBeta[1:10, ]
#'
#'
#' ## Using parameter coef and var.coef
#' saeHBbetaDeff <- BetaDeff(
#'   formula,
#'   ni = "ni", deff_iw = "deff",
#'   var.coef = vc, iter.update = 10,
#'   coef = c, data = dat
#' )
#'
#' saeHBbetaDeff$Est # Small Area mean Estimates
#' saeHBbetaDeff$refVar # Random effect variance
#' saeHBbetaDeff$coefficient # coefficient
#' # Load Library 'coda' to execute the plot
#' # autocorr.plot(saeHBbetaDeff$plot[[3]])  # is used to generate ACF Plot
#' # plot(saeHBbetaDeff$plot[[3]])           # is used to generate Density and trace plot
#'
#' ## Do not use parameter coef and var.coef
#' saeHBbetaDeff <- BetaDeff(formula, ni = "ni", deff_iw = "deff", data = dataBeta)
#' }
BetaDeff <- function(
    formula,
    ni = NULL,
    deff_iw = NULL,
    iter.update = 3, iter.mcmc = 10000,
    coef, var.coef, thin = 3, burn.in = 2000,
    tau.u = 1,
    data,
    seed = 1,
    quiet = TRUE,
    plot = TRUE) {
  result <- list(Est = NA, refVar = NA, coefficient = NA, plot = NA)
  formuladata <- stats::model.frame(formula, data, na.action = NULL)

  if (any(is.na(formuladata[, -1]))) {
    stop("Auxiliary Variables contains NA values.")
  }
  auxVar <- as.matrix(formuladata[, -1])
  nvar <- ncol(auxVar) + 1

  if (!missing(var.coef)) {
    if (length(var.coef) != nvar) {
      stop("length of vector var.coef does not match the number of regression coefficients, the length must be ", nvar)
    }
    tau.b.value <- 1 / var.coef
  } else {
    tau.b.value <- 1 / rep(1, nvar)
  }

  if (!missing(coef)) {
    if (length(coef) != nvar) {
      stop("length of vector coef does not match the number of regression coefficients, the length must be ", nvar)
    }
    mu.b.value <- coef
  } else {
    mu.b.value <- rep(0, nvar)
  }

  if (iter.update < 3) {
    stop("the number of iteration updates at least 3 times")
  }


  if (quiet) {
    cli::cli_progress_bar(
      total = iter.update,
      format = "Update {iter}/{iter.update} | {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} \n"
    )
    my_env <- environment()
    default_pb <- "none"
    update_progress <- function(iter, iter.update) {
      Sys.sleep(2 / 10000)
      cli::cli_progress_update(set = iter, .envir = my_env)
    }
  } else {
    default_pb <- "text"
    update_progress <- function(iter, iter.update) {
      cli::cli_h1("Update {iter}/{iter.update}")
    }
  }

  # rjags::load.module("lecuyer")
  if (!any(is.na(formuladata[, 1]))) {
    formuladata <- as.matrix(stats::na.omit(formuladata))

    if (any(formuladata[, 1] <= 0) || any(formuladata[, 1] >= 1)) {
      stop("response variable must be 0 < ", formula[2], " < 1")
    }

    x <- stats::model.matrix(formula, data = as.data.frame(formuladata))
    x <- as.matrix(x)
    n <- nrow(formuladata)
    mu.b <- mu.b.value
    tau.b <- tau.b.value
    tau.ub <- tau.ua <- a.var <- 1
    ni <- .get_variable(data, ni)
    deff_iw <- .get_variable(data, deff_iw)


    for (iter in 1:iter.update) {
      update_progress(iter, iter.update)

      dat <- list(
        "n" = n, "nvar" = nvar, "y" = formuladata[, 1], "x" = as.matrix(formuladata[, 2:nvar]), "mu.b" = mu.b,
        "tau.b" = tau.b, "tau.ua" = tau.ua, "tau.ub" = tau.ub,
        "ni" = ni, "deff_iw" = deff_iw
      )

      inits <- list(
        u = rep(0, n), b = mu.b, tau.u = tau.u,
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = seed
      )

      my_model <- "model{
					for (i in 1:n) {
							y[i] ~ dbeta(A[i],B[i])
							A[i] <- mu[i] * phi[i]
              B[i] <- (1-mu[i]) * phi[i]
              logit(mu[i]) <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)
							phi[i] <- (ni[i]/deff_iw[i] - 1)
					}

					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u
		  	}"

      jags.m <- rjags::jags.model(
        textConnection(my_model),
        data = dat,
        inits = inits, n.chains = 1, n.adapt = 500,
        quiet = quiet
      )

      params <- c("mu", "a.var", "b", "tau.u")
      samps <- rjags::coda.samples(
        jags.m, params,
        n.iter = iter.mcmc,
        thin = thin, seed = seed,
        progress.bar = default_pb
      )

      samps1 <- stats::window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps <- summary(samps1)
      a.var <- result_samps$statistics[1]
      beta <- result_samps$statistics[2:(nvar + 1), 1:2]
      for (i in 1:nvar) {
        mu.b[i] <- beta[i, 1]
        tau.b[i] <- 1 / (beta[i, 2]^2)
      }

      tau.ua <- result_samps$statistics[2 + nvar + n, 1]^2 / result_samps$statistics[2 + nvar + n, 2]^2
      tau.ub <- result_samps$statistics[2 + nvar + n, 1] / result_samps$statistics[2 + nvar + n, 2]^2
    }

    result_samps <- summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("b[", idx.b.varnames, "]"), pattern = " ", replacement = "")
    }
    result_mcmc <- samps1[, c(2:(nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var <- result_samps$statistics[1]

    beta <- result_samps$statistics[2:(nvar + 1), 1:2]
    rownames(beta) <- b.varnames

    mu <- result_samps$statistics[(nvar + 2):(1 + nvar + n), 1:2]

    Estimation <- data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2 + nvar + n), ])
    q_mu <- Quantiles[(nvar + 2):(nvar + 1 + n), ]
    q_beta <- (Quantiles[2:(nvar + 1), ])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation, q_mu)
    colnames(Estimation) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")
  } else {
    formuladata <- as.data.frame(formuladata)
    n <- nrow(formuladata)
    mu.b <- mu.b.value
    tau.b <- tau.b.value
    tau.ub <- tau.ua <- phi.aa <- phi.ab <- phi.ba <- phi.bb <- a.var <- 1
    formuladata$idx <- 1:n
    formuladata$ni <- .get_variable(data, ni)
    formuladata$deff <- .get_variable(data, deff_iw)
    data_sampled <- stats::na.omit(formuladata)

    if (any(data_sampled[, 1] <= 0) || any(data_sampled[, 1] >= 1)) {
      stop("response variable must be 0 < ", formula[2], " < 1")
    }

    data_nonsampled <- formuladata[-data_sampled$idx, ]
    r <- data_nonsampled$idx
    n1 <- nrow(data_sampled)
    n2 <- nrow(data_nonsampled)

    for (iter in 1:iter.update) {
      update_progress(iter, iter.update)

      dat <- list(
        "n1" = n1, "n2" = n2, "nvar" = nvar, "y_sampled" = data_sampled[, 1], "x_sampled" = as.matrix(data_sampled[, 2:nvar]),
        "x_nonsampled" = as.matrix(data_nonsampled[, 2:nvar]), "mu.b" = mu.b,
        "tau.b" = tau.b, "tau.ua" = tau.ua, "tau.ub" = tau.ub,
        "phi.aa" = phi.aa, "phi.ab" = phi.ab, "phi.ba" = phi.ba, "phi.bb" = phi.bb,
        "ni" = ni, "deff_iw" = deff_iw
      )

      inits <- list(
        u = rep(0, n1), v = rep(0, n2),
        b = mu.b, tau.u = tau.u,
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = seed
      )

      my_model <- "model {
					for (i in 1:n1) {
							y_sampled[i] ~ dbeta(A[i],B[i])
							A[i] <- mu[i] * phi[i]
              B[i] <- (1-mu[i]) * phi[i]
              logit(mu[i]) <- b[1] +  sum(b[2:nvar]*x_sampled[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)
							phi[i] <- (ni[i]/deff_iw[i] - 1)

					}

					for (j in 1:n2) {
				    	y.nonsampled[j] ~ dbeta(A.nonsampled[j],B.nonsampled[j])
							A.nonsampled[j] <- mu.nonsampled[j] * phi.nonsampled[j]
              B.nonsampled[j] <- (1-mu.nonsampled[j]) * phi.nonsampled[j]
              logit(mu.nonsampled[j]) <- mu.b[1] +sum(mu.b[2:nvar]*x_nonsampled[j,])+v[j]
					    v[j]~dnorm(0,tau.u)
					    phi.nonsampled[j] ~ dgamma(phi.a,phi.b)

          }

					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}

					phi.a ~ dgamma(phi.aa,phi.ab)
					phi.b ~ dgamma(phi.ba,phi.bb)
					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}"

      jags.m <- rjags::jags.model(
        file = textConnection(my_model), data = dat,
        inits = inits, n.chains = 1, n.adapt = 500,
        quiet = quiet
      )

      params <- c("mu", "mu.nonsampled", "a.var", "b", "phi.a", "phi.b", "tau.u")

      samps <- rjags::coda.samples(
        jags.m, params,
        n.iter = iter.mcmc,
        thin = thin,
        progress.bar = default_pb
      )
      samps1 <- stats::window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps <- summary(samps1)
      a.var <- result_samps$statistics[1]
      beta <- result_samps$statistics[2:(nvar + 1), 1:2]

      mu.b <- beta[, 1]
      tau.b <- 1 / (beta[, 2]^2)

      phi.aa <- result_samps$statistics[2 + nvar + n, 1]^2 / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ab <- result_samps$statistics[2 + nvar + n, 1] / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ba <- result_samps$statistics[3 + nvar + n, 1]^2 / result_samps$statistics[3 + nvar + n, 2]^2
      phi.bb <- result_samps$statistics[3 + nvar + n, 1] / result_samps$statistics[3 + nvar + n, 2]^2
      tau.ua <- result_samps$statistics[4 + nvar + n, 1]^2 / result_samps$statistics[4 + nvar + n, 2]^2
      tau.ub <- result_samps$statistics[4 + nvar + n, 1] / result_samps$statistics[4 + nvar + n, 2]^2
    }

    result_samps <- summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i - 1)
      b.varnames[i] <- str_replace_all(paste("b[", idx.b.varnames, "]"), pattern = " ", replacement = "")
    }
    result_mcmc <- samps1[, c(2:(nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var <- result_samps$statistics[1]

    beta <- result_samps$statistics[2:(nvar + 1), 1:2]
    rownames(beta) <- b.varnames

    mu <- result_samps$statistics[(nvar + 2):(1 + nvar + n1), 1:2]
    mu.nonsampled <- result_samps$statistics[(2 + nvar + n1):(1 + nvar + n), 1:2]

    Estimation <- matrix(rep(0, n), n, 2)
    Estimation[r, ] <- mu.nonsampled
    Estimation[-r, ] <- mu
    Estimation <- as.data.frame(Estimation)

    Quantiles <- as.data.frame(result_samps$quantiles[1:(3 + nvar + n), ])
    q_beta <- (Quantiles[2:(nvar + 1), ])
    q_mu <- Quantiles[(nvar + 2):(nvar + 1 + n1), ]
    q_mu.nonsampled <- (Quantiles[(2 + nvar + n1):(1 + nvar + n), ])
    q_Estimation <- matrix(0, n, 5)
    for (i in 1:5) {
      q_Estimation[r, i] <- q_mu.nonsampled[, i]
      q_Estimation[-r, i] <- q_mu[, i]
    }
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation, q_Estimation)
    colnames(Estimation) <- c("MEAN", "SD", "2.5%", "25%", "50%", "75%", "97.5%")
  }
  result$Est <- Estimation
  result$refVar <- a.var
  result$coefficient <- beta
  result$result_mcmc <- result_mcmc
  result$dat <- dat
  # result$inits <- inits
  result$formula <- formula
  result$formuladata <- formuladata
  class(result) <- "saehb"

  tryCatch(
    {
      if (plot) {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))

        graphics::par(mar = c(2, 2, 2, 2))
        coda::autocorr.plot(result_mcmc, col = "brown2", lwd = 2)
        plot(result_mcmc, col = "brown2", lwd = 2)
      }

      cli::cli_h1('Coefficient')
      stats::printCoefmat(beta)
      return(result)
    },
    error = function(cond) {
      cli::cli_alert_danger(cond)

      cli::cli_h1('Coefficient')
      stats::printCoefmat(beta)
      return(result)
    }
  )
  return(result)
}



.get_variable <- function(data, variable) {
  if (length(variable) == nrow(data)) {
    return(variable)
  } else if (methods::is(variable, "character")) {
    if (variable %in% colnames(data)) {
      variable <- data[[variable]]
    } else {
      cli::cli_abort('variable "{variable}" is not found in the data')
    }
  } else if (methods::is(variable, "formula")) {
    # extract column name (class character) from formula
    variable <- data[[all.vars(variable)]]
  } else {
    cli::cli_abort('variable "{variable}" is not found in the data')
  }
  return(variable)
}




.get_variable <- function(data, variable) {
  if (length(variable) == nrow(data)) {
    return(variable)
  } else if (methods::is(variable, "character")) {
    if (variable %in% colnames(data)) {
      variable <- data[[variable]]
    } else {
      cli::cli_abort('variable "{variable}" is not found in the data')
    }
  } else if (methods::is(variable, "formula")) {
    # extract column name (class character) from formula
    variable <- data[[all.vars(variable)]]
  } else {
    cli::cli_abort('variable "{variable}" is not found in the data')
  }
  return(variable)
}
