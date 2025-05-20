# Calculate epistasis p-values and other descriptives for
# timed epistasis models.
timed_epistasis_pvalue_calc = function(fit) {
  epistasis_lik_fn = fit@call$minuslogl
  orig_args = rlang::call_args(fit@call)
  orig_args = orig_args[! names(orig_args) == 'minuslogl']
  
  withCallingHandlers(
    {
      # Activate null-epistasis mode
      environment(epistasis_lik_fn)$null_epistasis = TRUE 
      curr_args = c(list(minuslogl = epistasis_lik_fn,
                         fixed = list(ces_A_on_B = 1e3, ces_B_on_A = 1e3)), orig_args)
      if(! is.null(curr_args$control$ndeps)) {
        curr_args$control$ndeps = curr_args$control$ndeps[1:2] # Assuming ndeps doesn't vary
      }
      null_fit = do.call(bbmle::mle2, curr_args)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
      if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
        # a little dangerous to muffle, but so far these warnings are
        # quite rare and have been harmless
        invokeRestart("muffleWarning") 
      }
    }
  )
  
  
  v1_simple_estimate = bbmle::coef(null_fit)['ces_A0']
  v2_simple_estimate = bbmle::coef(null_fit)['ces_B0']
  param_start = setNames(list(v1_simple_estimate, v2_simple_estimate, v1_simple_estimate, v2_simple_estimate),
                         bbmle::parnames(epistasis_lik_fn))
  
  # Turn null-epistasis mode back off
  environment(epistasis_lik_fn)$null_epistasis = FALSE
  
  param_init_v1 = list(ces_A0 = v1_simple_estimate, ces_B0 = 1e3, 
                       ces_A_on_B = v1_simple_estimate, ces_B_on_A = 1e3)
  
  withCallingHandlers(
    {
      curr_args = c(list(minuslogl = epistasis_lik_fn, 
                         fixed = list(ces_A0 = v1_simple_estimate, ces_A_on_B = v1_simple_estimate)),
                    orig_args)
      curr_args$start = param_init_v1
      if(! is.null(curr_args$control$ndeps)) {
        curr_args$control$ndeps = curr_args$control$ndeps[1:2] # Assuming ndeps doesn't vary
      }
      
      refit_v1 = do.call(bbmle::mle2, curr_args)
    },
    warning = function(w) {
      if (conditionMessage(w) %ilike% "some parameters are on the boundary") {
        invokeRestart("muffleWarning")
      }
    }
  )
  v1_chisquared = as.numeric(-2 * (bbmle::logLik(refit_v1) - bbmle::logLik(fit)))
  p_v1 =  pchisq(v1_chisquared, df = 2, lower.tail = F)
  
  
  param_init_v2 = list(ces_A0 = 1e3, ces_B0 = v2_simple_estimate,
                       ces_A_on_B = 1e3, ces_B_on_A = v2_simple_estimate)
  
  withCallingHandlers(
    {
      curr_args = c(list(minuslogl = epistasis_lik_fn, 
                         fixed = list(ces_B0 = v2_simple_estimate, ces_B_on_A = v2_simple_estimate)),
                    orig_args)
      curr_args$start = param_init_v2
      if(! is.null(curr_args$control$ndeps)) {
        curr_args$control$ndeps = curr_args$control$ndeps[1:2]
      }
      refit_v2 = do.call(bbmle::mle2, curr_args)
      
    },
    warning = function(w) {
      if (conditionMessage(w) %ilike% "some parameters are on the boundary") {
        invokeRestart("muffleWarning")
      }
    }
  )
  v2_chisquared = as.numeric(-2 * (bbmle::logLik(refit_v2) - bbmle::logLik(fit)))
  p_v2 = pchisq(v2_chisquared, df = 2, lower.tail = F)
  
  chisquared = as.numeric(-2 * (bbmle::logLik(null_fit) - bbmle::logLik(fit)))
  p_epistasis = pchisq(chisquared, df = 2, lower.tail = F)
  
  # Collect v1, v2 rates in the same sample order
  with_just_1 = environment(epistasis_lik_fn)$with_just_1
  with_just_2 = environment(epistasis_lik_fn)$with_just_2
  with_both = environment(epistasis_lik_fn)$with_both
  with_neither = environment(epistasis_lik_fn)$with_neither
  v1_rates_samples_with = c(with_just_1[[1]], with_both[[1]])
  v1_rates_samples_without = c(with_just_2[[1]], with_neither[[1]])
  v1_rates = c(v1_rates_samples_with, v1_rates_samples_without)
  v2_rates_samples_with = c(with_just_2[[2]], with_both[[2]])
  v2_rates_samples_without = c(with_just_1[[2]], with_neither[[2]])
  v2_rates = c(v2_rates_samples_with, v2_rates_samples_without)
  
  # Variables named for parallelism with pairwise_epistasis_lik().
  params = bbmle::coef(fit)
  A = params[1] * v1_rates
  B = params[2] * v2_rates
  A_on_B = params[3] * v1_rates
  B_on_A = params[4] * v2_rates
  
  # See pairwise_epistasis_lik()
  p_wt = exp(-1 * (A+B))
  p_A = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
  p_B = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))
  p_AB = 1 - p_wt - p_A - p_B
  expected_nAB_epistasis = sum(p_AB)
  
  # Under no-epistasis model, p(AB) = p(A) * p(B)
  A = v1_simple_estimate * v1_rates
  B = v2_simple_estimate * v2_rates
  expected_nAB_null = sum((1 - exp(-1 * A)) * (1 - exp(-1 * B)))
  
  
  return(list(p_A_change = p_v1, p_B_change = p_v2, p_epistasis = p_epistasis, 
              expected_nAB_epistasis = expected_nAB_epistasis, expected_nAB_null = expected_nAB_null,
              AB_epistatic_ratio = expected_nAB_epistasis / expected_nAB_null,
              ces_A_null = v1_simple_estimate, ces_B_null = v2_simple_estimate))
}
