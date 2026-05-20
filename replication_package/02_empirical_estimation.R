# =====================================================================
# 02_empirical_estimation.R
# ---------------------------------------------------------------------
# Apply the cross-sectionally augmented distributed-lag family to a
# real panel dataset stored in an Excel file, and report the long-run
# coefficients for the three estimators used in the empirical section:
#
#     CS-DL    : linear cross-sectionally augmented distributed lag
#     CS-NDL   : asymmetric (nonlinear) CS-DL with +/- partial sums
#     FCS-NDL  : Fourier-augmented CS-NDL (adds smooth structural change)
#
# WHAT THIS DOES
#   For each estimator it runs unit-by-unit OLS, averages the unit
#   coefficients (Mean Group), and reports the long-run elasticity,
#   the mean-group standard error, and a Wald test of long-run
#   symmetry for the asymmetric specifications.
#
# DATA FORMAT
#   The Excel file must contain four columns:
#
#        id        time       y            x
#     ---------------------------------------------
#     (unit)     (year)    (dependent)  (regressor)
#
#     id   : identifier for each cross-sectional unit (e.g., country)
#     time : time period (e.g., year)
#     y    : dependent variable (e.g., ln CO2 emissions per capita)
#     x    : explanatory variable  (e.g., ln energy use per capita)
#
# HOW TO RUN
#   1. Install R (https://cran.r-project.org/) and optionally RStudio.
#   2. Put your Excel file in the same folder as this script.
#   3. Edit the four lines under "USER SETTINGS" below.
#   4. Press "Source" in RStudio, OR run:
#         source("02_empirical_estimation.R")
#   5. Results print to the console and are written to
#      "empirical_results.csv".
#
# REQUIREMENTS
#   Base R plus one package (readxl), installed automatically below.
# =====================================================================


# ---------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cloud.r-project.org")
}
library(readxl)


# ---------------------------------------------------------------------
# 1. USER SETTINGS  --  edit these to match your data
# ---------------------------------------------------------------------
excel_file <- "data.xlsx"   # name of your Excel file
sheet_name <- 1             # sheet number (1) or name ("Sheet1")
n_fourier  <- 1             # number of Fourier frequencies (1 is typical)

# Column names inside your Excel file (change if yours differ)
col_id   <- "id"
col_time <- "time"
col_y    <- "y"
col_x    <- "x"


# ---------------------------------------------------------------------
# 2. Read the data
# ---------------------------------------------------------------------
raw <- as.data.frame(read_excel(excel_file, sheet = sheet_name))

dat <- data.frame(
  id   = raw[[col_id]],
  time = raw[[col_time]],
  y    = as.numeric(raw[[col_y]]),
  x    = as.numeric(raw[[col_x]])
)
dat <- dat[complete.cases(dat), ]
dat <- dat[order(dat$id, dat$time), ]

cat(sprintf("Loaded %d observations: %d units, time from %s to %s\n",
            nrow(dat), length(unique(dat$id)),
            as.character(min(dat$time)), as.character(max(dat$time))))


# ---------------------------------------------------------------------
# 3. Common inputs: cross-sectional averages and Fourier terms
# ---------------------------------------------------------------------
units    <- unique(dat$id)
N_units  <- length(units)
all_time <- sort(unique(dat$time))
T_len    <- length(all_time)

# Cross-sectional averages of y and x (CCE augmentation)
csa <- aggregate(cbind(y, x) ~ time, data = dat, FUN = mean)
names(csa) <- c("time", "y_bar", "x_bar")

# Fourier terms approximating smooth structural change
k       <- 1:n_fourier
tt      <- seq_along(all_time)
sin_mat <- sapply(k, function(kk) sin(2 * pi * kk * tt / T_len))
cos_mat <- sapply(k, function(kk) cos(2 * pi * kk * tt / T_len))
fourier_df <- data.frame(time = all_time,
                         sinF = rowSums(as.matrix(sin_mat)),
                         cosF = rowSums(as.matrix(cos_mat)))


# ---------------------------------------------------------------------
# 4. Estimation engine
#    model = "linear"     -> CS-DL    (one long-run coefficient)
#    model = "asymmetric" -> CS-NDL   (positive / negative)
#    model = "fourier"    -> FCS-NDL  (asymmetric + Fourier terms)
# ---------------------------------------------------------------------
mean_group <- function(b) {
  b <- b[is.finite(b)]
  est  <- mean(b)
  se   <- sd(b) / sqrt(length(b))
  c(estimate = est, std_error = se, t_value = est / se, n = length(b))
}

estimate_family <- function(model) {

  pos_hat <- numeric(N_units)
  neg_hat <- numeric(N_units)

  for (i in seq_len(N_units)) {

    di <- dat[dat$id == units[i], ]
    di <- merge(di, csa,        by = "time")
    di <- merge(di, fourier_df, by = "time")
    di <- di[order(di$time), ]

    if (model == "linear") {
      # CS-DL: y on x (level) + cross-sectional averages
      fit <- tryCatch(lm(y ~ x + y_bar + x_bar, data = di),
                      error = function(e) NULL)
      pos_hat[i] <- if (is.null(fit)) NA else coef(fit)["x"]
      neg_hat[i] <- NA  # no asymmetry in the linear model

    } else {
      # Asymmetric partial sums of x
      dx      <- c(di$x[1], diff(di$x))
      di$xpos <- cumsum(pmax(dx, 0))
      di$xneg <- cumsum(pmin(dx, 0))

      if (model == "asymmetric") {
        # CS-NDL: y on +/- partial sums + cross-sectional averages
        fit <- tryCatch(lm(y ~ xpos + xneg + y_bar + x_bar, data = di),
                        error = function(e) NULL)
      } else if (model == "fourier") {
        # FCS-NDL: CS-NDL plus Fourier terms
        fit <- tryCatch(
          lm(y ~ xpos + xneg + y_bar + x_bar + sinF + cosF, data = di),
          error = function(e) NULL)
      }
      pos_hat[i] <- if (is.null(fit)) NA else coef(fit)["xpos"]
      neg_hat[i] <- if (is.null(fit)) NA else coef(fit)["xneg"]
    }
  }

  out <- list(pos = mean_group(pos_hat))
  if (model != "linear") {
    out$neg <- mean_group(neg_hat)
    # Wald test of long-run symmetry: H0 theta_pos = theta_neg
    d  <- pos_hat - neg_hat
    d  <- d[is.finite(d)]
    md <- mean(d); sed <- sd(d) / sqrt(length(d))
    w  <- (md / sed)^2                  # ~ chi-squared(1) under H0
    out$wald <- c(stat = w, p_value = 1 - pchisq(w, df = 1))
  }
  out
}


# ---------------------------------------------------------------------
# 5. Run the three estimators
# ---------------------------------------------------------------------
res_csdl   <- estimate_family("linear")
res_csndl  <- estimate_family("asymmetric")
res_fcsndl <- estimate_family("fourier")


# ---------------------------------------------------------------------
# 6. Assemble and report
# ---------------------------------------------------------------------
fmt <- function(v) ifelse(is.na(v), "---", sprintf("%.4f", v))

results <- data.frame(
  statistic = c("LR(+)  estimate", "LR(+)  std error",
                "LR(-)  estimate", "LR(-)  std error",
                "Wald (symmetry) stat", "Wald p-value"),
  CS_DL = c(fmt(res_csdl$pos["estimate"]), fmt(res_csdl$pos["std_error"]),
            "---", "---", "---", "---"),
  CS_NDL = c(fmt(res_csndl$pos["estimate"]),  fmt(res_csndl$pos["std_error"]),
             fmt(res_csndl$neg["estimate"]),  fmt(res_csndl$neg["std_error"]),
             fmt(res_csndl$wald["stat"]),     fmt(res_csndl$wald["p_value"])),
  FCS_NDL = c(fmt(res_fcsndl$pos["estimate"]), fmt(res_fcsndl$pos["std_error"]),
              fmt(res_fcsndl$neg["estimate"]), fmt(res_fcsndl$neg["std_error"]),
              fmt(res_fcsndl$wald["stat"]),    fmt(res_fcsndl$wald["p_value"])),
  stringsAsFactors = FALSE
)

cat("\n============= EMPIRICAL ESTIMATION RESULTS =============\n")
cat("Long-run elasticities by estimator (Mean Group)\n\n")
print(results, row.names = FALSE)

cat("\nNotes:\n")
cat(" LR(+) and LR(-) are the long-run elasticities with respect to\n")
cat(" positive and negative shocks to x. CS-DL is linear (one coefficient).\n")
cat(" The Wald test has H0: LR(+) = LR(-); a p-value below 0.05 indicates\n")
cat(" significant long-run asymmetry.\n")

write.csv(results, "empirical_results.csv", row.names = FALSE)
cat("\nResults saved to empirical_results.csv\n")


# =====================================================================
# REFERENCE OUTPUT (paper application)
# ---------------------------------------------------------------------
# In the manuscript, this estimator family is applied to CO2 emissions
# per capita (y) and energy use per capita (x) for 37 OECD economies,
# 1990-2022. The long-run elasticities reported there are:
#
#                       CS-DL      CS-NDL     FCS-NDL
#   LR(+) Energy        0.308      0.283      0.122
#   LR(-) Energy        ---        0.314      0.295
#   Wald p-value        ---        [0.746]    [0.078]
#
# The positive long-run elasticity of CO2 with respect to energy use is
# positive and significant across specifications; the Wald test does not
# reject long-run symmetry at the 5% level. The example data.xlsx shipped
# with this repository is a small illustrative panel, so its numbers will
# differ from the paper; replace data.xlsx with the full panel built from
# the source in DATA_SOURCE.md to reproduce the paper's values.
# =====================================================================
