this_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(file.path(this_dir, "sim_mult_wave3_plus.R")) # assumes run_sim is in same directory as sim_mult_wave2.R
# and pwd is set to that directory
run_sim3 = function(Nsim=1000,nevents,looks_at,rho_target,...){

# alpha_one <- 0.025            # one-sided alpha used inside the sim
# HR0 <- 0.7                    # NI margin (VE > 0.3)
# hr_true <- 0.15               # true HR in your data-generating model (VE = 85%)



# ---- RUN THE SIMS ----
out <- vector("list", Nsim)
seeds <- sample.int(.Machine$integer.max, Nsim)
for (b in seq_len(Nsim)) {
  if (b %% 50 == 0) message("Completed ", b, " of ", Nsim)
  #print(paste0("Completed simulation ", b, " out of ", Nsim))
  #out[[b]] <- sim_mult_wave3(seedy = seeds[b],nevents=nevents,looks_at=looks_at,...)
  out[[b]] <- sim_mult_wave3_plus(seedy = seeds[b],nevents=nevents,looks_at=looks_at,rho_target = rho_target,...)
}

# ---- FILTER VALID RUNS ---- 
# should be all of them
ok = vapply(out, function(x) is.list(x) && isFALSE(x$error), logical(1))
n_valid = sum(ok)
if (n_valid == 0) stop("All runs errored; inspect function inputs/installation.")

# ---- POWER ESTIMATE ----
rej = vapply(out[ok], `[[`, logical(1), "reject_one_sided")
power_hat = mean(rej)

# ---- Monte Carlo precision (95% CI for power) ----
mc_se = sqrt(power_hat * (1 - power_hat) / n_valid)
power_ci = c(power_hat - 1.96 * mc_se, power_hat + 1.96 * mc_se)

# ---- Early stop rate and average look (diagnostics) ----
stopped = vapply(out[ok], function(x) isTRUE(x$stopped_early), logical(1))
early_stop_rate = mean(stopped)

stop_look = vapply(out[ok], function(x) x$stop_look, integer(1))
avg_stop_look = mean(stop_look, na.rm = TRUE)  # ~1,2,3 for looks at 10,15,20

TrialTime = vapply(out[ok], function(x) x$time_for_cases, double(1))
NumVaccinated = vapply(out[ok], function(x) x$total_vaccinated, double(1))
# ---- Summary ----
list(
  n_runs_total = Nsim,
  n_valid = n_valid,
  power = power_hat,
  power_95ci = power_ci,
  early_stop_rate = early_stop_rate,
  avg_stop_look = avg_stop_look,
  look_counts = table(factor(stop_look, levels = 1:3, labels = c("10ev","15ev","20ev"))),
  meanVacc = median(NumVaccinated),
  medTime = median(TrialTime)
)
}