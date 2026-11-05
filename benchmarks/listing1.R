# Copyright (c) 2026 Xiang Chen. Part of the spDBL R package.
# SPDX-License-Identifier: MIT

# Listing 1 of the MLOSS manuscript, character for character below this
# header. Run it with listing1_run.R, which records its output and run time.

library(spDBL)
data("dt_calibration"); d <- dt_calibration
em  <- emulator_learn(d$training$input, d$training$responses, d$Nx, d$Ny,
                      covariance = "sigma2R", episode_window = 2,
                      nsam = 20, seed = 1)               # closed form, no MCMC
cal <- calibration_prepare(em, d$lower, d$upper)        # design GP
calibration_validate(cal, d$validation$input, d$validation$responses)$field
fit <- calibration_fit(cal, d$field$z, d$field$U, chains = 4, iter = 2000,
                       cores = 4, init = "random", seed = 1)  # NUTS per draw
summary(fit)                                            # eta, R-hat, ESS
fields <- calibration_predict(fit)                      # y, u and z draws
