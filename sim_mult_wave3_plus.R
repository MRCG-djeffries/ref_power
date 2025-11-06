library(survival)
library(coxphf)

sim_mult_wave3_plus = function(
    nevents = 20,
    looks_at = c(10, 15, 20),
    HR0 = 0.7,
    waves = 1/rep(3,3),
    wave_attack = c(0.06, 0.06, 0.06),
    hr_true = 0.15,
    seroprev = 0.10,
    ## original per-subject frailty variance
    frailty_var = 0.5,
    ## OPTIONAL cluster frailty controls (default OFF)
    cluster_size = 5,
    rho_target = 0.3,
    cluster_frailty_var = 0,      # if >0, enables cluster frailty
    mapping = "second",
    alpha_final = 0.025,          # one-sided alpha
    seedy = 7563685,
    max_enrolled = 1000
) {
  set.seed(seedy)
  stopifnot(tail(looks_at, 1) == nevents)
  
  ## If cluster frailty effectively OFF, run original
  if (is.null(rho_target) && (is.null(cluster_frailty_var) || cluster_frailty_var == 0)) {
    ## ---- ORIGINAL BEHAVIOR (unchanged) ----
    info_frac = looks_at / nevents
    z_alpha   = qnorm(1 - alpha_final)
    z_bound   = z_alpha / sqrt(info_frac)
    
    enrolled = 0
    total_cases = 0
    cases_placebo = 0
    cases_vaccine = 0
    time_all = numeric(0)
    status_all = numeric(0)
    group_all = numeric(0)
    
    L = length(looks_at)
    z_seen        = rep(NA_real_, L)
    ve_lower_seen = rep(NA_real_, L)
    hr_seen       = rep(NA_real_, L)
    
    while (total_cases < nevents & enrolled < max_enrolled) {
      enrolled = enrolled + 1
      group  = rbinom(1, 1, 0.5)
      immune = rbinom(1, 1, seroprev) == 1
      if (frailty_var > 0) {
        frail = rgamma(1, shape = 1/frailty_var, scale = frailty_var)
      } else {
        frail = 1
      }
      t_ind = 0
      event_occured = FALSE
      for (i in seq_along(waves)) {
        lambda0 = -log(1 - wave_attack[i]) / waves[i]
        lambda  = lambda0 * frail
        if (group == 1) lambda = lambda * hr_true
        if (immune == 1) lambda = 0
        t_wave = rexp(1, rate = pmax(lambda, .Machine$double.eps))
        if (t_wave <= waves[i]) { t_ind = t_ind + t_wave; event_occured = TRUE; break }
        t_ind = t_ind + waves[i]
      }
      time_all   = c(time_all, t_ind)
      status_all = c(status_all, as.numeric(event_occured))
      group_all  = c(group_all, group)
      if (event_occured) {
        total_cases = total_cases + 1
        if (group == 0) cases_placebo = cases_placebo + 1 else cases_vaccine = cases_vaccine + 1
      }
      
      if (total_cases %in% looks_at) {
        k = match(total_cases, looks_at)
        events_by_arm = tapply(status_all, group_all, sum)
        if (any(events_by_arm == 0)) {
          df = data.frame(time_all = time_all, status_all = status_all, group_all = group_all)
          fit_look = tryCatch(coxphf(Surv(time_all, status_all) ~ group_all, data = df, firth = TRUE),
                              error = function(e) NULL)
          if (is.null(fit_look)) return(list(error = TRUE))
        } else {
          fit_look = tryCatch(coxph(Surv(time_all, status_all) ~ group_all), error = function(e) NULL)
          if (is.null(fit_look)) return(list(error = TRUE))
        }
        loghr_hat_look = as.numeric(coef(fit_look)[1])
        se_loghr_look  = sqrt(as.numeric(vcov(fit_look)[1,1]))
        z_stat_look    = (loghr_hat_look - log(HR0)) / se_loghr_look
        U_hr_one_look  = exp(loghr_hat_look + z_alpha * se_loghr_look)
        ve_lower_look  = 1 - U_hr_one_look
        hr_hat_look    = exp(loghr_hat_look)
        
        z_seen[k]        = z_stat_look
        ve_lower_seen[k] = ve_lower_look
        hr_seen[k]       = hr_hat_look
        
        if (k < length(looks_at) && z_stat_look < -z_bound[k]) {
          return(list(
            reject_one_sided = TRUE,
            stopped_early    = TRUE,
            stop_look        = k,
            z_boundary       = z_bound[k],
            z_stat           = z_stat_look,
            ve_hat           = 1 - hr_hat_look,
            ve_lower         = ve_lower_look,
            hr_hat           = hr_hat_look,
            z_seen           = z_seen,
            ve_lower_seen    = ve_lower_seen,
            hr_seen          = hr_seen,
            total_cases      = total_cases,
            cases_placebo    = cases_placebo,
            cases_vaccine    = cases_vaccine,
            total_enrolled   = enrolled,
            time_for_cases   = max(time_all[status_all == 1]),
            total_vaccinated = sum(group_all),
            error            = FALSE
          ))
        }
      }
    }
    
    if (total_cases < nevents) return(list(error = TRUE))
    
    k_last = length(looks_at)
    if (!is.na(z_seen[k_last])) {
      z_stat   = z_seen[k_last]
      hr_hat   = hr_seen[k_last]
      ve_lower = ve_lower_seen[k_last]
      ve_hat   = 1 - hr_hat
    } else {
      events_by_arm = tapply(status_all, group_all, sum)
      if (any(events_by_arm == 0)) {
        df = data.frame(time_all = time_all, status_all = status_all, group_all = group_all)
        fit = tryCatch(coxphf(Surv(time_all, status_all) ~ group_all, data = df, firth = TRUE),
                       error = function(e) NULL)
        if (is.null(fit)) return(list(error=TRUE))
      } else {
        fit = tryCatch(coxph(Surv(time_all, status_all) ~ group_all), error = function(e) NULL)
        if (is.null(fit)) return(list(error=TRUE))
      }
      loghr_hat = as.numeric(coef(fit)[1])
      se_loghr  = sqrt(as.numeric(vcov(fit)[1,1]))
      hr_hat    = exp(loghr_hat)
      z_stat    = (loghr_hat - log(HR0)) / se_loghr
      U_hr_one  = exp(loghr_hat + z_alpha * se_loghr)
      ve_hat    = 1 - hr_hat
      ve_lower  = 1 - U_hr_one
      
      z_seen[k_last]        = z_stat
      ve_lower_seen[k_last] = ve_lower
      hr_seen[k_last]       = hr_hat
    }
    reject_one_sided = (z_stat < -z_bound[k_last])
    
    return(list(
      reject_one_sided = reject_one_sided,
      stopped_early    = FALSE,
      stop_look        = k_last,
      z_boundary       = z_bound[k_last],
      z_stat = z_stat,
      ve_hat = ve_hat,
      ve_lower = ve_lower,
      hr_hat = hr_hat,
      z_seen = z_seen,
      ve_lower_seen = ve_lower_seen,
      hr_seen = hr_seen,
      total_cases = total_cases,
      cases_placebo = cases_placebo,
      cases_vaccine = cases_vaccine,
      total_enrolled = enrolled,
      time_for_cases = max(time_all[status_all == 1]),
      total_vaccinated = sum(group_all),
      error = FALSE)
    )
    ## ---- END ORIGINAL BEHAVIOR ----
  }
  
  ## ===== Cluster frailty path (ONLY when requested) =====
  mapping = match.arg(mapping, choices = c("engineering","first","second","exact"))
  info_frac = looks_at / nevents
  z_alpha   = qnorm(1 - alpha_final)
  z_bound   = z_alpha / sqrt(info_frac)
  
  Lambda_total = sum(-log(1 - wave_attack))
  
  theta_from_rho = function(Lambda, rho, method) {
    stopifnot(rho > 0, rho < 1, Lambda > 0)
    if (method == "first")       return(rho / Lambda) # crapola
    if (method == "second")      return(rho * (1 - Lambda) / Lambda)
    if (method == "engineering") return(rho / (Lambda * (1 - rho))) #crapola
    f = function(theta){
      a = (1 + theta * Lambda)^(-1/theta)
      b = (1 + 2 * theta * Lambda)^(-1/theta)
      (b - a^2) / (a * (1 - a)) - rho
    }
    uniroot(f, lower = 1e-9, upper = 100)$root
  }
  
  theta_cluster = if (!is.null(cluster_frailty_var) && cluster_frailty_var > 0) {
    cluster_frailty_var
  } else {
    theta_from_rho(Lambda_total, rho_target, mapping)
  }
  theta_ind = max(0, frailty_var)
  
  enrolled = 0L; total_cases = 0L
  cases_placebo = 0L; cases_vaccine = 0L
  time_all = numeric(0); status_all = numeric(0); group_all = numeric(0); cluster_all = integer(0)
  
  L = length(looks_at)
  z_seen = rep(NA_real_, L); ve_lower_seen = rep(NA_real_, L); hr_seen = rep(NA_real_, L)
  
  current_cluster_id = 0L; current_cluster_f = 1.0; within_left = 0L
  draw_gamma = function(theta) if (theta > 0) rgamma(1, shape = 1/theta, scale = theta) else 1.0
  
  while (total_cases < nevents & enrolled < max_enrolled) {
    
    if (within_left == 0L) {
      current_cluster_id = current_cluster_id + 1L
      current_cluster_f  = draw_gamma(theta_cluster)
      within_left        = cluster_size
    }
    
    enrolled   = enrolled + 1L
    within_left = within_left - 1L
    
    group  = rbinom(1, 1, 0.5)
    immune = rbinom(1, 1, seroprev) == 1
    
    indiv_f = draw_gamma(theta_ind)
    frail_total = current_cluster_f * indiv_f
    
    t_ind = 0; event_occured = FALSE
    for (i in seq_along(waves)) {
      lambda0 = -log(1 - wave_attack[i]) / waves[i]
      lambda  = lambda0 * frail_total
      if (group == 1) lambda = lambda * hr_true
      if (immune == 1) lambda = 0
      t_wave = rexp(1, rate = pmax(lambda, .Machine$double.eps))
      if (t_wave <= waves[i]) { t_ind = t_ind + t_wave; event_occured = TRUE; break }
      t_ind = t_ind + waves[i]
    }
    
    time_all   = c(time_all, t_ind)
    status_all = c(status_all, as.numeric(event_occured))
    group_all  = c(group_all, group)
    cluster_all = c(cluster_all, current_cluster_id)
    
    if (event_occured) {
      total_cases = total_cases + 1L
      if (group == 0) cases_placebo = cases_placebo + 1L else cases_vaccine = cases_vaccine + 1L
    }
    
    if (total_cases %in% looks_at) {
      k = match(total_cases, looks_at)
      events_by_arm = tapply(status_all, group_all, sum)
      
      ## Use robust SEs by cluster; if it fails (e.g., separation), Firth fallback
      df = data.frame(time_all, status_all, group_all, cluster_all)
      fit_look = tryCatch(
        coxph(Surv(time_all, status_all) ~ group_all + cluster(cluster_all), data = df,robust = TRUE),
        error = function(e) NULL)
      
      if (is.null(fit_look)) {
        df2 = data.frame(time_all, status_all, group_all)
        fit_look = tryCatch(coxphf(Surv(time_all, status_all) ~ group_all, data = df2, firth = TRUE,robust = TRUE),
                             error = function(e) NULL)
        if (is.null(fit_look)) return(list(error = TRUE))
        loghr_hat_look = as.numeric(coef(fit_look)[1])
        se_loghr_look  = sqrt(as.numeric(vcov(fit_look)[1,1]))
      } else {
        loghr_hat_look = as.numeric(coef(fit_look)["group_all"])
        se_loghr_look  = sqrt(vcov(fit_look)["group_all","group_all"])
      }
      
      z_stat_look   = (loghr_hat_look - log(HR0)) / se_loghr_look
      U_hr_one_look = exp(loghr_hat_look + z_alpha * se_loghr_look)
      ve_lower_look = 1 - U_hr_one_look
      hr_hat_look   = exp(loghr_hat_look)
      
      z_seen[k]        = z_stat_look
      ve_lower_seen[k] = ve_lower_look
      hr_seen[k]       = hr_hat_look
      
      if (k < length(looks_at) && z_stat_look < -z_bound[k]) {
        return(list(
          reject_one_sided = TRUE,
          stopped_early    = TRUE,
          stop_look        = k,
          z_boundary       = z_bound[k],
          z_stat           = z_stat_look,
          ve_hat           = 1 - hr_hat_look,
          ve_lower         = ve_lower_look,
          hr_hat           = hr_hat_look,
          z_seen           = z_seen,
          ve_lower_seen    = ve_lower_seen,
          hr_seen          = hr_seen,
          total_cases      = total_cases,
          cases_placebo    = cases_placebo,
          cases_vaccine    = cases_vaccine,
          total_enrolled   = enrolled,
          clusters_enrolled= current_cluster_id,
          time_for_cases   = if (any(status_all==1)) max(time_all[status_all == 1]) else NA_real_,
          total_vaccinated = sum(group_all),
          Lambda_total     = Lambda_total,
          theta_cluster_used = theta_cluster,
          theta_ind_used     = theta_ind,
          error            = FALSE
        ))
      }
    }
  }
  
  if (total_cases < nevents) return(list(error = TRUE))
  
  k_last = length(looks_at)
  if (!is.na(z_seen[k_last])) {
    z_stat   = z_seen[k_last]
    hr_hat   = hr_seen[k_last]
    ve_lower = ve_lower_seen[k_last]
    ve_hat   = 1 - hr_hat
  } else {
    df = data.frame(time_all, status_all, group_all, cluster_all)
    fit = tryCatch(
      coxph(Surv(time_all, status_all) ~ group_all + cluster(cluster_all), data = df),
      error = function(e) NULL)
    if (is.null(fit)) {
      df2 = data.frame(time_all, status_all, group_all)
      fit = tryCatch(coxphf(Surv(time_all, status_all) ~ group_all, data = df2, firth = TRUE),
                      error = function(e) NULL)
      if (is.null(fit)) return(list(error=TRUE))
      loghr_hat = as.numeric(coef(fit)[1])
      se_loghr  = sqrt(as.numeric(vcov(fit)[1,1]))
    } else {
      loghr_hat = as.numeric(coef(fit)["group_all"])
      se_loghr  = sqrt(vcov(fit)["group_all","group_all"])
    }
    hr_hat   = exp(loghr_hat)
    z_stat   = (loghr_hat - log(HR0)) / se_loghr
    U_hr_one = exp(loghr_hat + z_alpha * se_loghr)
    ve_hat   = 1 - hr_hat
    ve_lower = 1 - U_hr_one
    
    z_seen[k_last]        = z_stat
    ve_lower_seen[k_last] = ve_lower
    hr_seen[k_last]       = hr_hat
  }
  reject_one_sided = (z_stat < -z_bound[k_last])
  
  list(
    reject_one_sided   = reject_one_sided,
    stopped_early      = FALSE,
    stop_look          = k_last,
    z_boundary         = z_bound[k_last],
    z_stat             = z_stat,
    ve_hat             = 1 - hr_hat,
    ve_lower           = ve_lower,
    hr_hat             = hr_hat,
    z_seen             = z_seen,
    ve_lower_seen      = ve_lower_seen,
    hr_seen            = hr_seen,
    total_cases        = total_cases,
    cases_placebo      = cases_placebo,
    cases_vaccine      = cases_vaccine,
    total_enrolled     = enrolled,
    clusters_enrolled  = if (exists("current_cluster_id")) current_cluster_id else NA_integer_,
    time_for_cases     = if (any(status_all==1)) max(time_all[status_all == 1]) else NA_real_,
    total_vaccinated   = sum(group_all),
    Lambda_total       = sum(-log(1 - wave_attack)),
    theta_cluster_used = if (exists("theta_cluster")) theta_cluster else 0,
    theta_ind_used     = frailty_var,
    error = FALSE
  )
}
