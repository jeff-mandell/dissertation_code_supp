#' timed_pairwise_epistasis_lik
#'
#' For a pair of variants (or two groups of variants), creates a likelihood function for a
#' model of pairwise epistasis with a "strong mutation, weak selection" assumption and an estimate from MutationTimeR of which variant occurred after (ie. in the context of) the other. 
#' 
#' The arguments to this function (except v1_first and v2_first) are automatically supplied by \code{ces_epistasis()} and \code{ces_gene_epistasis()}. The user needs to provide v1_first and v2_first through the lik_args of ces_gene_epistasis. 
#' 
#' @param with_just_1 two-item list of  baseline rates in v1/v2 for tumors with mutation in just the first variant(s)
#' @param with_just_2 two-item list of  baseline rates in v1/v2 for tumors with mutation in just the second variant(s)
#' @param with_both two-item list of  baseline rates in v1/v2 for tumors with mutation in both
#' @param with_neither two-item list of  baseline rates in v1/v2 for tumors with mutation in neither two-item list of baseline rates in v1/v2 for tumors with mutation in neither
#' @param v1_first Samples that have a confirmed mutation ordering of v1 followed by v2 (must be subset of samples in with_both)
#' @param v2_first Samples that have a confirmed mutation ordering of v2 followed by v1 (must be subset of samples in with_both)
#' @export
#' @return A likelihood function
timed_pairwise_epistasis_lik  <- function(with_just_1, with_just_2, with_both, with_neither, 
                                          v1_first, v2_first, null_epistasis = FALSE) {
  
  if(length(intersect(v1_first, v2_first)) != 0) {
    stop("Sample(s) appear in both v1_first and v2_first.")
  }
  if(length(v1_first) != uniqueN(v1_first)) {
    stop('A sample appears more than once in v1_first.')
  }
  if(length(v2_first) != uniqueN(v2_first)) {
    stop('A sample appears more than once in v2_first.')
  }
  
  samples_not_in_both = setdiff(c(v1_first, v2_first), names(with_both[[1]]))
  if(length(samples_not_in_both) > 0) {
    if(length(samples_not_in_both) > 4) {
      additional_missing = length(samples_not_in_both) - 3
      samples_not_in_both = c(samples_not_in_both[1:3], paste0('(and ', additional_missing, ' more)'))
    } else {
          stop('There are samples listed in v1_first or v2_first that do not have both variants:\n',
         paste0(samples_not_in_both, collapse = ', '), '.')
    }
  }
  
  unspecified <- setdiff(names(with_both[[1]]), c(v1_first, v2_first))
  

  # Split up the with_both rate lists based on ordering
  with_both_1b2 = lapply(with_both, "[", v1_first)
  with_both_2b1 = lapply(with_both, "[", v2_first)
  with_both_unspecified = lapply(with_both, "[", unspecified)
  
  # Synthetic counts are a possibility (not currently implementing)
  # if(length(v1_first) == 0 || length(v2_first) == 0) {
  #   with_both_1b2[[1]] = c(with_both_1b2[[1]], c(synA = 1))
  #   with_both_1b2[[2]] = c(with_both_1b2[[2]], c(synA = 1e-5))
  #   with_both_2b1[[1]] = c(with_both_2b1[[1]], c(synB = 1e-5))
  #   with_both_2b1[[2]] = c(with_both_2b1[[2]], c(synB = 1))
  # }
  
  fn = function(par) {
    
    # sometimes the pars end up as NaNs or NAs, possibly because of inappropriate optimization techniques
    if(! all(is.finite(par))) {
      return(1e200)
    }
    
    # two points of discontinuity we need to account for
    if((par[3] == par[1] + par[2]) |
       (par[4] == par[1] + par[2])){return(1e200)}
    
    sum_log_lik <- 0
    
    if(! is.null(with_neither)) {
      # log(P{wt}) = -(A + B)
      A = par[1] * with_neither[[1]]
      B = par[2] * with_neither[[2]]
      ll = -1 * (A + B)
      sum_log_lik = sum_log_lik + sum(ll)  
    }
    
    
    if(! is.null(with_just_1)) {
      A = par[1] * with_just_1[[1]]
      B = par[2] * with_just_1[[2]]
      B_on_A = B
      if(! null_epistasis) {
        B_on_A = par[4] * with_just_1[[2]]
      }
      
      lik  = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
      sum_log_lik = sum_log_lik + sum(log(lik))   
    }
    
    
    if(! is.null(with_just_2)) {
      A = par[1] * with_just_2[[1]]
      B = par[2] * with_just_2[[2]]
      if(null_epistasis) {
        A_on_B = A
      } else {
        A_on_B = par[3] * with_just_2[[1]]
      }
      
      lik = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))
      sum_log_lik = sum_log_lik + sum(log(lik))
    }
    
    
    if(! is.null(with_both_1b2)) {
      A = par[1] * with_both_1b2[[1]]
      B = par[2] * with_both_1b2[[2]]

      if(null_epistasis) {
        B_on_A = B
      } else {
        B_on_A = par[4] * with_both_1b2[[2]]  
      }
      
      lik  = (A / (A + B)) * (1 - exp(-1 * (A + B)) - (((A + B) / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))))
      sum_log_lik = sum_log_lik + sum(log(lik))
    }
    
    
    if(! is.null(with_both_2b1)) {
      A = par[1] * with_both_2b1[[1]]
      B = par[2] * with_both_2b1[[2]]
      if(null_epistasis) {
        A_on_B = A
      } else {
        A_on_B = par[3] * with_both_2b1[[1]]
      }
      
      lik  = (B / (A + B)) * (1 - exp(-1 * (A + B)) - (((A + B) / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))))
      sum_log_lik = sum_log_lik + sum(log(lik))
    }
    
    
    if(! is.null(with_both_unspecified)) {
      A = par[1] * with_both_unspecified[[1]]
      B = par[2] * with_both_unspecified[[2]]
      if(null_epistasis) {
        A_on_B = A
        B_on_A = B
      } else {
        A_on_B = par[3] * with_both_unspecified[[1]]
        B_on_A = par[4] * with_both_unspecified[[2]]
      }
      
      # p_wt = exp(-1 * (A+B))
      # p_A = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
      # p_B = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))   
      # p_AB_unspecified = 1 - p_wt - p_A - p_B
      # sum_log_lik = sum_log_lik + sum(log(p_AB_unspecified))
      
      lik = (B / (A + B)) * (1 - exp(-1 * (A + B)) - (((A + B) / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B))))) +
        (A / (A + B)) * (1 - exp(-1 * (A + B)) - (((A + B) / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))))
      sum_log_lik = sum_log_lik + sum(log(lik))
    }
    
    
    # in case it tried all the max at once.
    if(!is.finite(sum_log_lik)){
      return(1e200)
    }
    return(-sum_log_lik)
  }
  
  # Set default values for all parameters, which ces_variant will use to set starting values of optimization
  formals(fn)[["par"]] = 1000:1003
  
  # Optimization tool, bbmle::mle, requires that vector of parameters to optimize have named elements
  bbmle::parnames(fn) = c("ces_A0", "ces_B0", "ces_A_on_B", "ces_B_on_A")
  return(fn)
}







