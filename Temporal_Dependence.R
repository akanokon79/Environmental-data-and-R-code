# Temporal Dependence Diagnostic
# Regularized Quantile Regression (QR / RQR / ARQR)
# Purpose: Examine within-country temporal dependence in residuals
# Diagnostics: Lag-1 Residual Autocorrelation (ACF1) & Ljung-Box Test (Lag 3)
# Quantiles: 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95
# Target Countries: Benin, Ghana, Nigeria, Senegal, Togo


rm(list = ls())

# 1. Environment & Dependencies 
if (!requireNamespace("quantreg", quietly = TRUE)) {
  install.packages("quantreg")
}
library(quantreg)

# Reproducibility seed
set.seed(2026)

# 2. Data Ingestion & Preprocessing 
data_path <- "data/evdata.csv"
dat <- read.csv(data_path, stringsAsFactors = FALSE)

# Handle missingness and specify country factor
dat <- na.omit(dat)
names(dat)[1] <- "Country"
dat$Country <- as.factor(dat$Country)

# Chronological sorting within country panel
dat <- dat[order(dat$Country, dat$Year), , drop = FALSE]

# Standardize variable naming conventions
names(dat) <- gsub("\\(", ".", names(dat))
names(dat) <- gsub("\\)", ".", names(dat))

# Set reference level (Benin)
if ("Benin" %in% levels(dat$Country)) {
  dat$Country <- factor(dat$Country, levels = c("Benin", setdiff(levels(dat$Country), "Benin")))
} else {
  warning("Benin was not found. The first country will be used as the reference.")
}

# Define response and continuous predictor specifications
y <- as.numeric(dat$LNC02)
continuous_names <- c("LNGDPp", "GDPg", "FDI", "LN.EU.", "LN.EC.", "AE", "REC", "UPG", "TNR")

# Construct design matrix components
X_cont <- as.matrix(dat[, continuous_names])
storage.mode(X_cont) <- "numeric"

country_dummies <- model.matrix(~ Country, data = dat)[, -1, drop = FALSE]
X_raw <- as.matrix(cbind(X_cont, country_dummies))
storage.mode(X_raw) <- "numeric"

# 3. Model Configuration & Selected Hyperparameters 
taus <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
models <- c("QR", "RQR", "ARQR")

# Adaptive parameter specifications
gamma <- 1
eps <- 1e-6

# Dimension tracking
p_cont <- length(continuous_names)
p_total <- ncol(X_raw)
p_fe <- p_total - p_cont
param_names <- c("Intercept", colnames(X_raw))

# Selected lambda tuning values from revision analysis
selected_lambda <- data.frame(
  Tau = c(0.05, 0.05, 0.10, 0.10, 0.25, 0.25, 0.50, 0.50, 0.75, 0.75, 0.90, 0.90, 0.95, 0.95),
  Model = rep(c("RQR", "ARQR"), 7),
  Selected_Lambda = c(0.25, 10, 0.25, 10, 2, 25, 5, 75, 5, 75, 5, 25, 5, 25)
)

# Lambda grid used in the revised analysis
lambda_grid <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50,
                 1, 2, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500, 750, 1000)

# 4. Core Statistical Functions 

quantile_loss <- function(y, yhat, tau) {
  u <- y - yhat
  mean(u * (tau - (u < 0)))
}

# Objective function for Penalized Quantile Regression
penalized_qr_objective <- function(beta, X, y, tau, lambda, penalty_weights) {
  fitted <- as.vector(X %*% beta)
  residuals <- y - fitted
  check_loss <- sum(residuals * (tau - (residuals < 0)))
  ridge_penalty <- (lambda / 2) * sum(penalty_weights * beta^2)
  return(check_loss + ridge_penalty)
}

# Unpenalized Quantile Regression (rq.fit.fnb)
fit_QR <- function(X, y, tau) {
  Xp <- cbind(Intercept = 1, X)
  fit <- rq.fit.fnb(Xp, y, tau = tau)
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  names(beta) <- colnames(Xp)
  yhat <- as.vector(Xp %*% beta)
  list(beta = beta, yhat = yhat)
}

# Calculate ARQR adaptive weights: w_j = (|beta_j^QR| + eps)^(-gamma)
calculate_adaptive_weights <- function(X, y, tau, gamma = 1, eps = 1e-6) {
  Xp <- cbind(Intercept = 1, X)
  qr_fit <- rq.fit.fnb(Xp, y, tau = tau)
  beta_qr <- qr_fit$coefficients
  beta_qr[is.na(beta_qr)] <- 0
  
  beta_cont <- beta_qr[2:(p_cont + 1)]
  raw_weights <- (abs(beta_cont) + eps)^(-gamma)
  geometric_mean <- exp(mean(log(raw_weights)))
  adaptive_weights <- raw_weights / geometric_mean
  
  list(
    beta_qr = beta_qr,
    raw_weights = raw_weights,
    adaptive_weights = adaptive_weights
  )
}

# Penalized QR Estimation Wrapper via Optimization
fit_penalized <- function(X, y, tau, model, lambda, gamma = 1, eps = 1e-6, adaptive_weights_external = NULL) {
  Xp <- cbind(Intercept = 1, X)
  qr_fit <- rq.fit.fnb(Xp, y, tau = tau)
  beta_qr <- qr_fit$coefficients
  beta_qr[is.na(beta_qr)] <- 0
  
  if (model == "RQR") {
    continuous_weights <- rep(1, p_cont)
  } else if (model == "ARQR") {
    if (is.null(adaptive_weights_external)) {
      weight_info <- calculate_adaptive_weights(X = X, y = y, tau = tau, gamma = gamma, eps = eps)
      continuous_weights <- weight_info$adaptive_weights
    } else {
      continuous_weights <- adaptive_weights_external
    }
  }
  
  # Penalty profile: Intercept = 0, Predictors = continuous_weights, Fixed Effects = 0
  penalty_weights <- c(0, continuous_weights, rep(0, p_fe))
  
  # Multi-start initializations for numerical stability
  starting_values <- list(
    beta_qr,
    beta_qr * 0.90,
    beta_qr * 0.75,
    beta_qr * 0.50,
    beta_qr * 0.25
  )
  
  best_value <- Inf
  best_beta <- beta_qr
  
  for (start in starting_values) {
    opt <- try(
      optim(
        par = start,
        fn = penalized_qr_objective,
        X = Xp,
        y = y,
        tau = tau,
        lambda = lambda,
        penalty_weights = penalty_weights,
        method = "Nelder-Mead",
        control = list(maxit = 50000, reltol = 1e-10)
      ),
      silent = TRUE
    )
    
    if (!inherits(opt, "try-error") && all(is.finite(opt$par))) {
      if (opt$value < best_value) {
        best_value <- opt$value
        best_beta <- opt$par
      }
    }
  }
  
  beta <- best_beta
  names(beta) <- colnames(Xp)
  yhat <- as.vector(Xp %*% beta)
  
  list(
    beta = beta,
    yhat = yhat,
    penalty_weights = penalty_weights,
    objective = best_value
  )
}

# Stratified Sample Split Utility
stratified_sample <- function(country, proportion) {
  selected <- c()
  for (cc in levels(country)) {
    idx <- which(country == cc)
    n_select <- floor(proportion * length(idx))
    if (n_select < 1) n_select <- 1
    selected <- c(selected, sample(idx, size = n_select, replace = FALSE))
  }
  return(selected)
}

# 5. Data Partitioning & Normalization 

# Stratified 80% Training sample selection
train_full_idx <- stratified_sample(dat$Country, 0.80)
X_train_full_raw <- X_raw[train_full_idx, , drop = FALSE]
y_train_full <- y[train_full_idx]

# Continuous predictor normalization
X_train_full <- X_train_full_raw
train_scaled <- scale(X_train_full[, 1:p_cont, drop = FALSE])
train_center <- attr(train_scaled, "scaled:center")
train_scale  <- attr(train_scaled, "scaled:scale")

X_train_full[, 1:p_cont] <- train_scaled

# 6. Model Estimation & Temporal Diagnostics 

temporal_results <- data.frame()
coefficient_results <- data.frame()

for (tau in taus) {
  for (model in models) {
    cat("\n\n")
    cat("TEMPORAL DIAGNOSTIC:", model, "| TAU =", tau, "\n")
    cat("\n")
    
    if (model == "QR") {
      lambda_final <- 0
      final_weights <- NULL
      fit_final <- fit_QR(X = X_train_full, y = y_train_full, tau = tau)
    } else {
      lambda_final <- selected_lambda[
        selected_lambda$Tau == tau & selected_lambda$Model == model,
        "Selected_Lambda"
      ]
      
      if (model == "ARQR") {
        final_weight_info <- calculate_adaptive_weights(
          X = X_train_full, y = y_train_full, tau = tau, gamma = gamma, eps = eps
        )
        final_weights <- final_weight_info$adaptive_weights
      } else {
        final_weights <- NULL
      }
      
      fit_final <- fit_penalized(
        X = X_train_full,
        y = y_train_full,
        tau = tau,
        model = model,
        lambda = lambda_final,
        gamma = gamma,
        eps = eps,
        adaptive_weights_external = final_weights
      )
    }
    
    # Store estimated coefficients
    coefficient_results <- rbind(
      coefficient_results,
      data.frame(
        Tau = tau,
        Model = model,
        Lambda = lambda_final,
        Parameter = names(fit_final$beta),
        Estimate = as.numeric(fit_final$beta)
      )
    )
    
    # Extract and order residuals by country and year
    residuals_final <- y_train_full - fit_final$yhat
    
    diagnostic_data <- data.frame(
      Country = dat$Country[train_full_idx],
      Year = dat$Year[train_full_idx],
      Residual = residuals_final
    )
    diagnostic_data <- diagnostic_data[
      order(diagnostic_data$Country, diagnostic_data$Year), , drop = FALSE
    ]
    
    # Within-country temporal autocorrelation & independence tests
    for (cc in levels(dat$Country)) {
      country_data <- diagnostic_data[diagnostic_data$Country == cc, , drop = FALSE]
      residual_country <- country_data$Residual
      n_country <- length(residual_country)
      
      # ACF(1) calculation
      if (n_country >= 2) {
        acf_result <- acf(residual_country, lag.max = 1, plot = FALSE, na.action = na.pass)
        acf1 <- as.numeric(acf_result$acf[2])
      } else {
        acf1 <- NA_real_
      }
      
      # Ljung-Box test through Lag 3
      if (n_country > 3) {
        lb_result <- Box.test(residual_country, lag = 3, type = "Ljung-Box")
        lb_stat   <- as.numeric(lb_result$statistic)
        lb_pvalue <- as.numeric(lb_result$p.value)
      } else {
        lb_stat   <- NA_real_
        lb_pvalue <- NA_real_
      }
      
      # Flag dependence based on significance threshold
      if (!is.na(lb_pvalue)) {
        dependence_flag <- if (lb_pvalue < 0.05) "Evidence of temporal dependence" else "No significant evidence"
      } else {
        dependence_flag <- NA_character_
      }
      
      temporal_results <- rbind(
        temporal_results,
        data.frame(
          Model = model,
          Tau = tau,
          Country = cc,
          N = n_country,
          ACF1 = acf1,
          Ljung_Box_Lag = 3,
          Ljung_Box_Statistic = lb_stat,
          Ljung_Box_P_Value = lb_pvalue,
          Dependence_Flag = dependence_flag
        )
      )
    }
  }
}

# 7. Output Formatting & Tabulation 

temporal_results$ACF1 <- round(temporal_results$ACF1, 4)
temporal_results$Ljung_Box_Statistic <- round(temporal_results$Ljung_Box_Statistic, 4)
temporal_results$Ljung_Box_P_Value <- signif(temporal_results$Ljung_Box_P_Value, 5)

cat("\n\n\n")
cat("TEMPORAL DEPENDENCE DIAGNOSTIC RESULTS\n")
cat("\n\n")
print(temporal_results, row.names = FALSE)

# Aggregate summary of autocorrelation violations
temporal_summary <- aggregate(
  Ljung_Box_P_Value ~ Model + Tau,
  data = temporal_results,
  FUN = function(x) sum(x < 0.05, na.rm = TRUE)
)
names(temporal_summary)[names(temporal_summary) == "Ljung_Box_P_Value"] <- "Countries_with_p_less_than_0.05"

cat("\n\n\n")
cat("SUMMARY OF TEMPORAL DEPENDENCE\n")
cat("\n\n")
print(temporal_summary, row.names = FALSE)

# 8. Export Results 

write.csv(temporal_results, "Temporal_Dependence_QR_RQR_ARQR_7Quantiles_FINAL.csv", row.names = FALSE)
write.csv(coefficient_results, "Temporal_Diagnostic_Model_Coefficients_QR_RQR_ARQR_7Quantiles.csv", row.names = FALSE)
write.csv(temporal_summary, "Temporal_Dependence_Summary_QR_RQR_ARQR_7Quantiles.csv", row.names = FALSE)

cat("\n\n\n")
cat("TEMPORAL DEPENDENCE DIAGNOSTIC COMPLETED\n")
cat("\n")
cat("Models: QR, RQR, ARQR\n")
cat("Quantiles: 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95\n")
cat("Countries: Benin, Ghana, Nigeria, Senegal, Togo\n")
cat("Diagnostic: ACF(1) + Ljung-Box test through lag 3\n")
cat("Outputs Saved:\n")
cat(" - Temporal_Dependence_QR_RQR_ARQR_7Quantiles_FINAL.csv\n")
cat(" - Temporal_Diagnostic_Model_Coefficients_QR_RQR_ARQR_7Quantiles.csv\n")
cat(" - Temporal_Dependence_Summary_QR_RQR_ARQR_7Quantiles.csv\n")
cat("\n")