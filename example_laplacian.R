rm(list=ls())
hablar::set_wd_to_script_path()

# laplacian prior
prior_density <- function(x){
  0.5*exp(-abs(x))
}
prior_cdf <- function(x){
  ifelse(x<0, 0.5*exp(x),
         1 -0.5*exp(-x))
}

# the encoding function is the cumulative of the prior (i.e., pnorm())
# the inverse-encoding (from internal to stimulus space) is given by quantile function
# to calculate likelihood function in stimulus space we need to aply change of variable and 
# take the jacobian into account (the derivative of the encoding function, which is the prob density function)
likelihood_s <- function(x, x_internal, noise_sd=0.05)
{
  dens <- rep(NA,length(x))
  for(i in 1:length(x)){
    jacobian <- abs(prior_density(x[i]))
    dens[i] <- dnorm(prior_cdf(x[i]), mean=x_internal, sd=noise_sd) * jacobian
  }
  return(dens)
}

# # sanity check
FUN <- function(x){ likelihood_s(x, x_internal=pnorm(0)) }
integrate(FUN, lower=0, upper= Inf)

# plot

# grid
x <- seq(-4,4, length.out=1e4)
prior_x <- prior_density(x)


# plot likelihood functions 
lik_0 <- likelihood_s(x, x_internal=prior_cdf(0))
lik_1 <- likelihood_s(x, x_internal=prior_cdf(1))
lik_2 <- likelihood_s(x, x_internal=prior_cdf(2))

# scale for plotting visibility
prior_x <- prior_x/max(prior_x) * max(lik_0)/2

palette <- viridis::inferno(4)

png("laplacian_prior.png", width=800, height=800)

par(mfrow=c(2,2))
plot(x, prior_x, col="dark grey", type="l", lwd=4, ylab = "density", ylim=c(0, max(lik_0)), main="example likelihoods")
lines(x, lik_0, col= palette[1], lwd=4); abline(v=0, col=palette[1])
lines(x, lik_1, col= palette[2], lwd=4); abline(v=1, col=palette[2])
lines(x, lik_2, col= palette[3], lwd=4); abline(v=2, col=palette[3])
legend('topright',"prior",lwd=4, col="dark grey",bty="n")


# calculate posterior mean for a bunch of values
# (setting x_internal to expected value)
x_i <- seq(-2,2, length.out=50)
post_mean <- rep(NA, 50)
post_sd <- rep(NA, 50)
for(i in 1:length(x_i)){
  lik_i <- likelihood_s(x, x_internal=pnorm(x_i[i]))
  post_i <-  (lik_i * prior_x) / sum(lik_i * prior_x)
  post_mean[i] <- sum(x*post_i)
  
  # estimate SD via sampling
  post_sample <- sample(x, size = 10000, prob = post_i, replace = TRUE)
  post_sd[i] <- sd(post_sample)
}
plot(x_i, post_mean, type="l",lwd=3, col="blue",xlab="true value",ylab="expected posterior mean")
abline(a=0,b=1, lty=2)

plot(x_i, post_mean-x_i, type="l",lwd=3,col="blue",xlab="true value",ylab="bias")
abline(h=0, lty=2)
abline(v=0, lty=2)

plot(x_i, post_sd, type="l",lwd=3,col="blue",xlab="true value",ylab="posterior SD")

dev.off()

