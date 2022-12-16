
# Code adapted from Douma & Weedo 2019 - Appendix S4
# Need to run 3_group_NOTIME.Rmd before

library(mvtnorm)

# function to combine the computation of the mean and the confidence limits
mean.quant <- function(x, probs) {
  a <- mean(x, na.rm = T)
  b <- t(quantile(x, probs = probs, na.rm = T))
  d <- c(a, b)
  names(d) <- c("mean", "lower", "upper")
  return(d)
}


####### GROW #######

means <- unlist(coef(grow2))
vc <- vcov(grow2)

N=1000

rnd <- rmvnorm(N, mean = means, sigma = vc)

modobj <- grow2 # calls best model fit object `modobj`

dm <- do.call(cbind, modobj$X) # is the design matrix, don't need the last column
dp <- modobj$Z
n <- nrow(diridatag)

# nao da pra pegar por nomes, estou pegando por posicao
quad_lp <- rnd[ ,1:2] %*% t(dm)[1:2, ]
sp_lp <- rnd[ ,3:4] %*% t(dm)[3:4, ] 
quadsp_lp <- rnd[ ,5:6] %*% t(dm)[5:6, ] 
res_lp <- matrix(0, ncol = n, nrow = N) 

precision_lp <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp) # precision

## convert linear predictors to mus
mu_output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))


for (i in 1:N) {
  denom <- (exp(quad_lp[i, ]) + exp(sp_lp[i, ])+ exp(quadsp_lp[i, ]) + 
              exp(rep(0, n)))
  mu_output[, "simQUAD", i] <- exp(quad_lp[i, ]) / denom
  mu_output[, "simSP", i] <- exp(sp_lp[i, ]) / denom
  mu_output[, "simQUADSP", i] <- exp(quadsp_lp[i, ]) / denom
  mu_output[, "simRES", i] <- 1 - (mu_output[, "simQUAD", i] + 
                                     mu_output[, "simSP", i] + 
                                     mu_output[, "simQUADSP", i])
}

# calculate alphas
alphas_QUAD <- mu_output[, "simQUAD", ] * t(exp(precision_lp))
alphas_SP <- mu_output[, "simSP", ] * t(exp(precision_lp))
alphas_QUADSP <- mu_output[, "simQUADSP", ] * t(exp(precision_lp))
alphas_RES <- mu_output[, "simRES", ] * t(exp(precision_lp))


# simulate from data 
output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))

for (i in 1:dim(alphas_RES)[2]) {
  # prediction
  test <- cbind(alphas_QUAD[, i], alphas_SP[, i], alphas_QUADSP[, i],
                alphas_RES[,i])
  output[, , i] <- rdirichlet(n, test)
}


# generate mean and 95% quantile limits
QUAD_quant <- apply(output[, "simQUAD", ], 1, mean.quant, probs = c(0.025, 0.975))
SP_quant <- apply(output[, "simSP", ], 1, mean.quant, probs = c(0.025, 0.975))
QUADSP_quant <- apply(output[, "simQUADSP", ], 1, mean.quant, probs = c(0.025, 0.975))
RES_quant <- apply(output[, "simRES", ], 1, mean.quant, probs = c(0.025, 0.975))

# add  confidence limits to data.frame
quant_QUAD <- data.frame(diridatag[, c("log.rich")], "term" = "quadrat", 
                         t(QUAD_quant))
quant_SP <- data.frame(diridatag[, c("log.rich")], "term" = "sp", 
                       t(SP_quant))
quant_QUADSP <- data.frame(diridatag[, c("log.rich")], "term" = "quadrat:sp", t(QUADSP_quant))
quant_RES <- data.frame(diridatag[, c("log.rich")], "term" = "Residual", t(RES_quant))
quantg <- rbind.data.frame(quant_QUAD, quant_SP, quant_QUADSP, quant_RES)

quantg$log.rich.o <- quantg$log.rich*sd(diridatag$log.rich.o) +
  mean(diridatag$log.rich.o)
quantg$rich.o <- exp(quantg$log.rich.o)


# newpredg <- newpredg %>% mutate(term = fct_relevel(term, "quadrat", "quadrat:sp", "sp",
#                                                    "Residual"))
# # ggplot(newpredg, aes(x=rich.o, y=pred)) +
#   geom_line() + facet_grid(~term) +
#   geom_point(data=predirig, aes(x=richness.rarefaction, y=VPC, col=term)) +
#   ggtitle("Mod log.rich") +
#   scale_x_log10() +
#   geom_smooth(data=quantg,aes(x=rich.o, y=mean), col="red", se=F)+
#   geom_smooth(data=quantg,aes(x=rich.o, y=lower), se=F)+
#   geom_smooth(data=quantg,aes(x=rich.o, y=upper), se=F)  +
#   geom_ribbon(data=quantg,aes(x=rich.o, ymin=lower, ymax=upper, y=mean), alpha=0.1, fill="blue",
#               size=0)




####### MORT #######

means <- unlist(coef(mort2))
vc <- vcov(mort2)

N=1000

rnd <- rmvnorm(N, mean = means, sigma = vc)

modobj <- mort2 # calls best model fit object `modobj`

dm <- do.call(cbind, modobj$X) # is the design matrix, don't need the last column
dp <- modobj$Z
n <- nrow(diridatam)

# nao da pra pegar por nomes, estou pegando por posicao
quad_lp <- rnd[ ,1:2] %*% t(dm)[1:2, ]
sp_lp <- rnd[ ,3:4] %*% t(dm)[3:4, ] 
quadsp_lp <- rnd[ ,5:6] %*% t(dm)[5:6, ] 
res_lp <- matrix(0, ncol = n, nrow = N) 

precision_lp <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp) # precision

## convert linear predictors to mus
mu_output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))


for (i in 1:N) {
  denom <- (exp(quad_lp[i, ]) + exp(sp_lp[i, ])+ exp(quadsp_lp[i, ]) + 
              exp(rep(0, n)))
  mu_output[, "simQUAD", i] <- exp(quad_lp[i, ]) / denom
  mu_output[, "simSP", i] <- exp(sp_lp[i, ]) / denom
  mu_output[, "simQUADSP", i] <- exp(quadsp_lp[i, ]) / denom
  mu_output[, "simRES", i] <- 1 - (mu_output[, "simQUAD", i] + 
                                     mu_output[, "simSP", i] + 
                                     mu_output[, "simQUADSP", i])
}

# calculate alphas
alphas_QUAD <- mu_output[, "simQUAD", ] * t(exp(precision_lp))
alphas_SP <- mu_output[, "simSP", ] * t(exp(precision_lp))
alphas_QUADSP <- mu_output[, "simQUADSP", ] * t(exp(precision_lp))
alphas_RES <- mu_output[, "simRES", ] * t(exp(precision_lp))


# simulate from data 
output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))

for (i in 1:dim(alphas_RES)[2]) {
  # prediction
  test <- cbind(alphas_QUAD[, i], alphas_SP[, i], alphas_QUADSP[, i],
                alphas_RES[,i])
  output[, , i] <- rdirichlet(n, test)
}


# generate mean and 95% quantile limits
QUAD_quant <- apply(output[, "simQUAD", ], 1, mean.quant, probs = c(0.025, 0.975))
SP_quant <- apply(output[, "simSP", ], 1, mean.quant, probs = c(0.025, 0.975))
QUADSP_quant <- apply(output[, "simQUADSP", ], 1, mean.quant, probs = c(0.025, 0.975))
RES_quant <- apply(output[, "simRES", ], 1, mean.quant, probs = c(0.025, 0.975))

# add  confidence limits to data.frame
quant_QUAD <- data.frame(diridatam[, c("log.rich")], "term" = "quadrat", 
                         t(QUAD_quant))
quant_SP <- data.frame(diridatam[, c("log.rich")], "term" = "sp", 
                       t(SP_quant))
quant_QUADSP <- data.frame(diridatam[, c("log.rich")], "term" = "quadrat:sp", t(QUADSP_quant))
quant_RES <- data.frame(diridatam[, c("log.rich")], "term" = "Residual", t(RES_quant))
quantm <- rbind.data.frame(quant_QUAD, quant_SP, quant_QUADSP, quant_RES)


quantm$log.rich.o <- quantm$log.rich*sd(diridatam$log.rich.o) +
  mean(diridatam$log.rich.o)
quantm$rich.o <- exp(quantm$log.rich.o)




# teste <- predirim %>% left_join(diridatam[,c(1,3,11)], by=c("fplot"))


# newpredm <- newpredm %>% mutate(term = fct_relevel(term, "quadrat", "quadrat:sp", "sp",
#                                                    "Residual"))
# ggplot(newpredm, aes(x=log.rich, y=pred)) +
#   geom_line() + facet_grid(~term) +
#   geom_point(data=teste, aes(x=log.rich, y=VPC, col=term)) +
#   ggtitle("Mod log.rich") +
#   geom_smooth(data=quantm,aes(x=log.rich, y=mean), col="red", se=F)+
#   geom_smooth(data=quantm,aes(x=log.rich, y=lower), se=F)+
#   geom_smooth(data=quantm,aes(x=log.rich, y=upper), se=F) 





####### REC #######

means <- unlist(coef(rec2))
vc <- vcov(rec2)

N=1000

rnd <- rmvnorm(N, mean = means, sigma = vc)

modobj <- rec2 # calls best model fit object `modobj`

dm <- do.call(cbind, modobj$X) # is the design matrix, don't need the last column
dp <- modobj$Z
n <- nrow(diridatar)

# nao da pra pegar por nomes, estou pegando por posicao
quad_lp <- rnd[ ,1:2] %*% t(dm)[1:2, ]
sp_lp <- rnd[ ,3:4] %*% t(dm)[3:4, ] 
quadsp_lp <- rnd[ ,5:6] %*% t(dm)[5:6, ] 
res_lp <- matrix(0, ncol = n, nrow = N) 

precision_lp <- rnd[, grep("gamma", colnames(rnd))] %*% t(dp) # precision

## convert linear predictors to mus
mu_output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))


for (i in 1:N) {
  denom <- (exp(quad_lp[i, ]) + exp(sp_lp[i, ])+ exp(quadsp_lp[i, ]) + 
              exp(rep(0, n)))
  mu_output[, "simQUAD", i] <- exp(quad_lp[i, ]) / denom
  mu_output[, "simSP", i] <- exp(sp_lp[i, ]) / denom
  mu_output[, "simQUADSP", i] <- exp(quadsp_lp[i, ]) / denom
  mu_output[, "simRES", i] <- 1 - (mu_output[, "simQUAD", i] + 
                                     mu_output[, "simSP", i] + 
                                     mu_output[, "simQUADSP", i])
}

# calculate alphas
alphas_QUAD <- mu_output[, "simQUAD", ] * t(exp(precision_lp))
alphas_SP <- mu_output[, "simSP", ] * t(exp(precision_lp))
alphas_QUADSP <- mu_output[, "simQUADSP", ] * t(exp(precision_lp))
alphas_RES <- mu_output[, "simRES", ] * t(exp(precision_lp))


# simulate from data 
output <- array(NA, dim = c(n, 4, N), dimnames = list(c(1:n), c("simQUAD", "simSP", "simQUADSP", "simRES"), c(1:N)))

for (i in 1:dim(alphas_RES)[2]) {
  # prediction
  test <- cbind(alphas_QUAD[, i], alphas_SP[, i], alphas_QUADSP[, i],
                alphas_RES[,i])
  output[, , i] <- rdirichlet(n, test)
}


# generate mean and 95% quantile limits
QUAD_quant <- apply(output[, "simQUAD", ], 1, mean.quant, probs = c(0.025, 0.975))
SP_quant <- apply(output[, "simSP", ], 1, mean.quant, probs = c(0.025, 0.975))
QUADSP_quant <- apply(output[, "simQUADSP", ], 1, mean.quant, probs = c(0.025, 0.975))
RES_quant <- apply(output[, "simRES", ], 1, mean.quant, probs = c(0.025, 0.975))

# add  confidence limits to data.frame
quant_QUAD <- data.frame(diridatar[, c("log.rich")], "term" = "quadrat", 
                         t(QUAD_quant))
quant_SP <- data.frame(diridatar[, c("log.rich")], "term" = "sp", 
                       t(SP_quant))
quant_QUADSP <- data.frame(diridatar[, c("log.rich")], "term" = "quadrat:sp", t(QUADSP_quant))
quant_RES <- data.frame(diridatar[, c("log.rich")], "term" = "Residual", t(RES_quant))
quantr <- rbind.data.frame(quant_QUAD, quant_SP, quant_QUADSP, quant_RES)

quantr$log.rich.o <- quantr$log.rich*sd(diridatar$log.rich.o) +
  mean(diridatar$log.rich.o)
quantr$rich.o <- exp(quantr$log.rich.o)



#teste <- predirir %>% left_join(diridatar[,c(1,11)], by=c("fplot"))


# newpredr <- newpredr %>% mutate(term = fct_relevel(term, "quadrat", "quadrat:sp", "sp",
#                                                    "Residual"))
# ggplot(newpredr, aes(x=log.rich, y=pred)) +
#   geom_line() + facet_grid(~term) +
#   geom_point(data=teste, aes(x=log.rich, y=VPC, col=term)) +
#   ggtitle("Mod log.rich") +
#   geom_smooth(data=quantr,aes(x=log.rich, y=mean), col="red", se=F)+
#   geom_smooth(data=quantr,aes(x=log.rich, y=lower), se=F)+
#   geom_smooth(data=quantr,aes(x=log.rich, y=upper), se=F) 
# 
# 


####### saving RESULTS ######

#save(quantg,quantm,quantr, file="results/prediction_intervals_dirichlet_exclude.Rdata")
