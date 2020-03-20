#' To simulate spatially coorelated cpue data over the sampling area
#' 
#' @param N The "true" number of fish
#' @param sim_locs A data.frame representation of the coordinate locations over which cpue data will be simulated 
#' @param variog_fit A vector containing estimates of the nugget, partial sill, and range parameter for the exponential variogram. The parameter estimates must be given in that order.
#' @param q_bar An estimate of the catchability coefficient
#' @param sd_q The estimated standard deviation of the catchability coefficient
#' @param n_sim The number of simulations
#' @param A The area of the lake (in XX)
#' @param n_max For local kriging -- he number of nearest observations that should be used for kriging simulation. see help(gstat)
#' @return Returns a data.frame containing the coordinate locations simulated cpue values
#' @export

cpue.sim <- function(N, sim_locs, variog_fit, q_bar, sd_q, n_sim, A, n_max=100){
  q <- rnorm(n_sim, q_hat, q_bar)
  cpue <- N/A*q
  gs <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=cpue, model=vgm(nugget=variog_fit[1], psill=variog_fit[2], range=variog_fit[3], model='Exp'), nmax=n_max)
  yy <- predict(gs, newdata=sim_locs, nsim=n_sim)
  return(yy)
}