# Regularized Quantile Regression for Assessing Heterogeneous Carbon Emissions
# Predictors in Selected West African Countries
# Models: Standard Quantile Regression (QR), Ridge QR (RQR), Adaptive Ridge QR (ARQR)
# Target Quantiles: 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95
# Inference: Within-country 3-year moving block bootstrap (B = 500)

rm(list = ls())

# 1. Environment & Dependencies 
if (!requireNamespace("quantreg", quietly = TRUE)) {
  stop("The 'quantreg' package is required. Install it via install.packages('quantreg').")
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
  warning("Benin not found in factor levels. Defaulting to first observed level as reference.")
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

# 3. Model Configuration & Hyperparameters 
taus <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
models <- c("QR", "RQR", "ARQR")

# Adaptive parameter specifications
gamma <- 1
eps <- 1e-6

# Penalty tuning grid
lambda_grid <- c(
  0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 
  1, 2, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500, 750, 1000
)

# Resampling settings
B <- 500

# Dimension tracking
p_cont <- length(continuous_names)
p_total <- ncol(X_raw)
p_fe <- p_total - p_cont
param_names <- c("Intercept", colnames(X_raw))

# 4. Core Statistical Functions 
# Check/Quantile loss evaluation
quantile_loss <- function(y, yhat, tau) {
  u <- y - yhat
  mean(u * (tau - (u < 0)))
}

# Objective function for Penalized Quantile Regression
# Q(beta) = sum rho_tau(y_i - x_i'beta) + (lambda / 2) * sum(w_j * beta_j^2)
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
  
  # Multi-start initializations for robust numerical optimization
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

# Stratified Sample Split Utility Function
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

# Initial Train (80%) and Test (20%) Split
train_full_idx <- stratified_sample(dat$Country, 0.80)
test_idx <- setdiff(seq_len(nrow(dat)), train_full_idx)

# Sub-split Estimation (75% of Train) and Validation (25% of Train)
country_train_full <- dat$Country[train_full_idx]
est_relative_idx <- stratified_sample(country_train_full, 0.75)
est_idx <- train_full_idx[est_relative_idx]
validation_idx <- setdiff(train_full_idx, est_idx)

# Data subsets extraction
X_est_raw <- X_raw[est_idx, , drop = FALSE]
y_est <- y[est_idx]

X_validation_raw <- X_raw[validation_idx, , drop = FALSE]
y_validation <- y[validation_idx]

X_train_full_raw <- X_raw[train_full_idx, , drop = FALSE]
y_train_full <- y[train_full_idx]

X_test_raw <- X_raw[test_idx, , drop = FALSE]
y_test <- y[test_idx]

# Scaling: Estimation & Validation split
X_est <- X_est_raw
X_validation <- X_validation_raw

est_scaled <- scale(X_est[, 1:p_cont, drop = FALSE])
est_center <- attr(est_scaled, "scaled:center")
est_scale  <- attr(est_scaled, "scaled:scale")

X_est[, 1:p_cont] <- est_scaled
X_validation[, 1:p_cont] <- scale(
  X_validation[, 1:p_cont, drop = FALSE],
  center = est_center,
  scale = est_scale
)

# 6. Hyperparameter Tuning (Lambda Selection) 

lambda_results <- data.frame()
selected_lambda <- data.frame()
master_coefficients <- data.frame()
master_performance <- data.frame()
master_weights <- data.frame()

for (tau in taus) {
  cat("\n\n\n")
  cat("LAMBDA TUNING: TAU =", tau, "\n")
  cat("\n")
  
  for (model in c("RQR", "ARQR")) {
    cat("\nTuning", model, "\n")
    
    if (model == "ARQR") {
      weight_info <- calculate_adaptive_weights(
        X = X_est, y = y_est, tau = tau, gamma = gamma, eps = eps
      )
      tuning_weights <- weight_info$adaptive_weights
    } else {
      tuning_weights <- NULL
    }
    
    validation_loss <- numeric(length(lambda_grid))
    
    for (k in seq_along(lambda_grid)) {
      lambda_value <- lambda_grid[k]
      
      fit_lambda <- fit_penalized(
        X = X_est,
        y = y_est,
        tau = tau,
        model = model,
        lambda = lambda_value,
        gamma = gamma,
        eps = eps,
        adaptive_weights_external = tuning_weights
      )
      
      Xp_validation <- cbind(Intercept = 1, X_validation)
      yhat_validation <- as.vector(Xp_validation %*% fit_lambda$beta)
      
      validation_loss[k] <- quantile_loss(y_validation, yhat_validation, tau)
      
      lambda_results <- rbind(
        lambda_results,
        data.frame(
          Tau = tau,
          Model = model,
          Lambda = lambda_value,
          Validation_Quantile_Loss = validation_loss[k]
        )
      )
      
      cat("Lambda =", lambda_value, "| Validation QL =", round(validation_loss[k], 6), "\n")
    }
    
    best_position <- which.min(validation_loss)
    best_lambda <- lambda_grid[best_position]
    best_validation_loss <- validation_loss[best_position]
    
    selected_lambda <- rbind(
      selected_lambda,
      data.frame(
        Tau = tau,
        Model = model,
        Selected_Lambda = best_lambda,
        Validation_Quantile_Loss = best_validation_loss
      )
    )
    
    cat("\nSELECTED:", model, "| lambda =", best_lambda, "\n")
  }
}

cat("\n\n\n")
cat("SELECTED LAMBDA VALUES\n")
cat("\n")
print(selected_lambda, row.names = FALSE)

# 7. Final Normalization & Model Fitting 

X_train_full <- X_train_full_raw
X_test <- X_test_raw

train_scaled <- scale(X_train_full[, 1:p_cont, drop = FALSE])
train_center <- attr(train_scaled, "scaled:center")
train_scale  <- attr(train_scaled, "scaled:scale")

X_train_full[, 1:p_cont] <- train_scaled
X_test[, 1:p_cont]       <- scale(X_test[, 1:p_cont, drop = FALSE], center = train_center, scale = train_scale)

for (tau in taus) {
  for (model in models) {
    cat("\n\n\n")
    cat("FINAL MODEL:", model, "| TAU =", tau, "\n")
    cat("\n")
    
    if (model == "QR") {
      fit_final <- fit_QR(X = X_train_full, y = y_train_full, tau = tau)
      lambda_final <- 0
      final_weights <- NULL
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
        
        master_weights <- rbind(
          master_weights,
          data.frame(
            Tau = tau,
            Model = model,
            Variable = continuous_names,
            QR_Estimate = round(final_weight_info$beta_qr[2:(p_cont + 1)], 6),
            Raw_Weight = round(final_weight_info$raw_weights, 6),
            Adaptive_Weight = round(final_weights, 6)
          )
        )
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
    
    # Prediction Generation
    yhat_train <- fit_final$yhat
    
    Xp_test <- cbind(Intercept = 1, X_test)
    yhat_test <- as.vector(Xp_test %*% fit_final$beta)
    
    # In-Sample and Out-of-Sample Evaluation
    train_mae   <- mean(abs(y_train_full - yhat_train))
    train_qloss <- quantile_loss(y_train_full, yhat_train, tau)
    test_mae    <- mean(abs(y_test - yhat_test))
    test_qloss  <- quantile_loss(y_test, yhat_test, tau)
    
    # 8. Within-Country Moving-Block Bootstrap Inference ----------------------
    block_length <- 3
    train_row_ids <- train_full_idx
    
    train_boot_meta <- data.frame(
      RowID = train_row_ids,
      Country = dat$Country[train_row_ids],
      Year = dat$Year[train_row_ids]
    )
    
    train_boot_meta <- train_boot_meta[
      order(train_boot_meta$Country, train_boot_meta$Year), , drop = FALSE
    ]
    
    # Moving-block bootstrap sampler function
    moving_block_sample <- function(n, block_length) {
      if (block_length > n) {
        stop("Block length cannot exceed the number of observations.")
      }
      possible_starts <- 1:(n - block_length + 1)
      n_blocks <- ceiling(n / block_length)
      starts <- sample(possible_starts, size = n_blocks, replace = TRUE)
      
      local_index <- unlist(lapply(starts, function(s) s:(s + block_length - 1)), use.names = FALSE)
      local_index[seq_len(n)]
    }
    
    country_groups_boot <- split(train_boot_meta$RowID, train_boot_meta$Country)
    boot_beta <- matrix(NA_real_, nrow = B, ncol = length(param_names))
    
    for (b in seq_len(B)) {
      boot_idx_original <- unlist(
        lapply(country_groups_boot, function(row_ids) {
          local_idx <- moving_block_sample(
            n = length(row_ids),
            block_length = min(block_length, length(row_ids))
          )
          row_ids[local_idx]
        }),
        use.names = FALSE
      )
      
      boot_positions <- match(boot_idx_original, train_full_idx)
      boot_positions <- boot_positions[!is.na(boot_positions)]
      
      X_b <- X_train_full[boot_positions, , drop = FALSE]
      y_b <- y_train_full[boot_positions]
      
      if (model == "ARQR") {
        boot_weight_info <- calculate_adaptive_weights(
          X = X_b, y = y_b, tau = tau, gamma = gamma, eps = eps
        )
        boot_weights <- boot_weight_info$adaptive_weights
      } else {
        boot_weights <- NULL
      }
      
      boot_fit <- try(
        if (model == "QR") {
          fit_QR(X = X_b, y = y_b, tau = tau)
        } else {
          fit_penalized(
            X = X_b,
            y = y_b,
            tau = tau,
            model = model,
            lambda = lambda_final,
            gamma = gamma,
            eps = eps,
            adaptive_weights_external = boot_weights
          )
        },
        silent = TRUE
      )
      
      if (!inherits(boot_fit, "try-error")) {
        if (length(boot_fit$beta) == length(param_names)) {
          boot_beta[b, ] <- boot_fit$beta
        }
      }
      
      if (b %% 50 == 0) {
        cat("Block bootstrap", b, "/", B, "\n")
      }
    }
    
    # Percentile-based 95% Confidence Intervals
    ci_lower <- apply(boot_beta, 2, function(z) {
      if (all(is.na(z))) return(NA_real_)
      quantile(z, probs = 0.025, na.rm = TRUE, names = FALSE)
    })
    
    ci_upper <- apply(boot_beta, 2, function(z) {
      if (all(is.na(z))) return(NA_real_)
      quantile(z, probs = 0.975, na.rm = TRUE, names = FALSE)
    })
    
    # Store coefficient outputs
    coefficient_table <- data.frame(
      Tau = tau,
      Model = model,
      Lambda = lambda_final,
      Parameter = param_names,
      Estimate = round(fit_final$beta, 4),
      CI_Lower_95 = round(ci_lower, 4),
      CI_Upper_95 = round(ci_upper, 4)
    )
    master_coefficients <- rbind(master_coefficients, coefficient_table)
    
    # Store evaluation metrics
    performance_table <- data.frame(
      Tau = tau,
      Model = model,
      Lambda = lambda_final,
      Train_MAE = round(train_mae, 4),
      Train_Quantile_Loss = round(train_qloss, 4),
      Test_MAE = round(test_mae, 4),
      Test_Quantile_Loss = round(test_qloss, 4)
    )
    master_performance <- rbind(master_performance, performance_table)
    
    cat("\nSelected lambda:", lambda_final, "\n")
    cat("Test MAE:", round(test_mae, 6), "\n")
    cat("Test Quantile Loss:", round(test_qloss, 6), "\n")
  }
}

# 9. Output Consolidation & Tabulation 
options(max.print = 100000)

cat("\n\n\n")
cat("FINAL COEFFICIENT ESTIMATES WITH 95% BOOTSTRAP CIs\n")
cat("\n\n")
print(master_coefficients, row.names = FALSE)
cat("\nTotal coefficient rows:", nrow(master_coefficients), "\n")
cat("Total coefficient columns:", ncol(master_coefficients), "\n")

cat("\n\n\n")
cat("QR: COMPLETE COEFFICIENT ESTIMATES\n")
cat("\n\n")
print(subset(master_coefficients, Model == "QR"), row.names = FALSE)

cat("\n\n\n")
cat("RQR: COMPLETE COEFFICIENT ESTIMATES\n")
cat("\n\n")
print(subset(master_coefficients, Model == "RQR"), row.names = FALSE)

cat("\n\n\n")
cat("ARQR: COMPLETE COEFFICIENT ESTIMATES\n")
cat("\n\n")
print(subset(master_coefficients, Model == "ARQR"), row.names = FALSE)

cat("\n\n\n")
cat("FINAL OUT-OF-SAMPLE PREDICTIVE PERFORMANCE\n")
cat("\n\n")
print(master_performance, row.names = FALSE)

cat("\n\n\n")
cat("FINAL ARQR ADAPTIVE WEIGHTS\n")
cat("\n\n")
print(master_weights, row.names = FALSE)

# Minimum observed test loss by quantile
lowest_observed_ql <- do.call(
  rbind,
  lapply(taus, function(tau_value) {
    tmp <- master_performance[master_performance$Tau == tau_value, , drop = FALSE]
    tmp[which.min(tmp$Test_Quantile_Loss), , drop = FALSE]
  })
)

cat("\n\n\n")
cat("LOWEST OBSERVED TEST QUANTILE LOSS BY QUANTILE\n")
cat("\n\n")
print(lowest_observed_ql, row.names = FALSE)

# Comparative metrics relative to standard QR
qr_reference <- master_performance[
  master_performance$Model == "QR",
  c("Tau", "Test_MAE", "Test_Quantile_Loss")
]

comparison_results <- merge(master_performance, qr_reference, by = "Tau", suffixes = c("", "_QR"))

comparison_results$MAE_Percent_Change_vs_QR <- 100 * 
  (comparison_results$Test_MAE_QR - comparison_results$Test_MAE) / comparison_results$Test_MAE_QR

comparison_results$QL_Percent_Change_vs_QR <- 100 * 
  (comparison_results$Test_Quantile_Loss_QR - comparison_results$Test_Quantile_Loss) / comparison_results$Test_Quantile_Loss_QR

comparison_results <- comparison_results[comparison_results$Model != "QR", ]

cat("\n\n\n")
cat("DESCRIPTIVE COMPARISON OF RQR AND ARQR WITH QR\n")
cat("\n\n")
print(
  comparison_results[, c(
    "Tau", "Model", "Lambda", "Test_MAE", 
    "Test_Quantile_Loss", "MAE_Percent_Change_vs_QR", "QL_Percent_Change_vs_QR"
  )],
  row.names = FALSE
)

# 10. Result Export -----------------------------------------------------------
write.csv(selected_lambda, "Selected_Lambda_QR_RQR_ARQR_FINAL.csv", row.names = FALSE)
write.csv(lambda_results, "Lambda_Validation_QR_RQR_ARQR_FINAL.csv", row.names = FALSE)
write.csv(master_coefficients, "Coefficients_QR_RQR_ARQR_FINAL.csv", row.names = FALSE)
write.csv(master_performance, "Performance_QR_RQR_ARQR_FINAL.csv", row.names = FALSE)
write.csv(master_weights, "Adaptive_Weights_ARQR_FINAL.csv", row.names = FALSE)
write.csv(comparison_results, "QR_RQR_ARQR_Percent_Comparison_FINAL.csv", row.names = FALSE)
write.csv(lowest_observed_ql, "Lowest_Observed_Test_QL_By_Quantile_FINAL.csv", row.names = FALSE)

cat("\n\n\n")
cat("RESULT FILES SAVED\n")
cat("\n")
cat("Selected_Lambda_QR_RQR_ARQR_FINAL.csv\n")
cat("Lambda_Validation_QR_RQR_ARQR_FINAL.csv\n")
cat("Coefficients_QR_RQR_ARQR_FINAL.csv\n")
cat("Performance_QR_RQR_ARQR_FINAL.csv\n")
cat("Adaptive_Weights_ARQR_FINAL.csv\n")
cat("QR_RQR_ARQR_Percent_Comparison_FINAL.csv\n")
cat("Lowest_Observed_Test_QL_By_Quantile_FINAL.csv\n")
cat("\n")

cat("\n\n\n")
cat("QR / RQR / ARQR ANALYSIS ACROSS SEVEN QUANTILES COMPLETED\n")
cat("\n")