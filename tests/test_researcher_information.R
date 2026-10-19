# The complete choice layer runs in an isolated environment without the DGP,
# config, menu, oracle, population functions, or the global workspace.
stopifnot(identical(parent.env(environment(researcher_choice)),baseenv()),
          !exists("config",envir=environment(researcher_choice),inherits=TRUE),
          !exists("oracle_loss",envir=environment(researcher_choice),inherits=TRUE))
e <- new.env(parent=baseenv())
sys.source("R/researcher_information.R",e)
sys.source("R/perceived_loss.R",e)
for (name in c("config","oracle_loss","conditional_mean","population_estimands","gamma","sigma","theta_a"))
  makeActiveBinding(name,function(value) stop("Hidden truth was accessed"),e)
p <- config$dgp
dat <- simulate_sample(p,config$seeds$test+200L)
obs <- dat[,c("Y","X","D","G")]
a <- e$researcher_choice(obs)
stopifnot(a$decision$ok, all(c(a$decision$perceived_loss_2,a$decision$perceived_loss_5)>=0))
# Extra truth columns are refused. A caller's evaluation truth cannot change choice.
stopifnot(inherits(try(e$researcher_choice(dat),silent=TRUE),"try-error"))
poison <- dat; poison$m <- NA_real_; poison$ps <- -1
poison$gamma <- 999; poison$oracle_loss <- -Inf; poison$theta_a <- 888; poison$sigma <- 777
stopifnot(identical(a,e$researcher_choice(poison[,c("Y","X","D","G")])) )
# Controlled deterministic augmented mean; coefficient and nested-model OVB identity.
obs$Y <- 1+.4*obs$X+.2*obs$X^2+.3*obs$G+.25*obs$D+.25*obs$D*obs$G+.17*obs$G*obs$X
b <- e$researcher_choice(obs)
expect_close(b$information$gamma_hat,.17)
expect_close(b$information$T_2,.5)
expect_close(b$information$T_augmented,.5)
expect_close(b$information$estimated_bias_5,b$information$T_5-b$information$T_augmented)
# Independent lm checks for noisy data, classical contrast variance and HC3.
obs <- dat[,c("Y","X","D","G")]
aug <- lm(Y~D+I(D*G)+X+I(X^2)+G+I(G*X),data=obs)
five <- lm(Y~D+I(D*G)+X+I(X^2)+G,data=obs)
two <- lm(Y~D+X+I(X^2),data=obs,subset=G==1)
expect_close(a$information$gamma_hat,coef(aug)["I(G * X)"])
expect_close(a$information$T_5,coef(five)["D"]+coef(five)["I(D * G)"])
expect_close(a$information$estimated_bias_5,a$information$T_5-(coef(aug)["D"]+coef(aug)["I(D * G)"]))
W <- model.matrix(five); contrast <- setNames(rep(0,ncol(W)),colnames(W))
contrast[c("D","I(D * G)")] <- 1
s2 <- sum(resid(aug)^2)/df.residual(aug)
expect_close(a$information$estimated_variance_5,drop(s2*t(contrast)%*%solve(crossprod(W),contrast)))
expect_close(a$information$estimated_variance_2,s2*solve(crossprod(model.matrix(two)))["D","D"])
expect_close(a$information$hc3_variance_2,fit_analysis(dat,2)$SE^2)
expect_close(a$information$hc3_variance_5,fit_analysis(dat,5)$SE^2)
# Changing treatment levels with terms in both models leaves the design rule
# invariant; dependence on outcome noise still remains through correlated scores.
shifted <- obs; shifted$Y <- shifted$Y+3*shifted$D-2*shifted$D*shifted$G
z <- e$researcher_choice(shifted)
expect_close(c(z$decision$perceived_loss_2,z$decision$perceived_loss_5),
             c(a$decision$perceived_loss_2,a$decision$perceived_loss_5))
fixture <- a$information; fixture$estimated_bias_5 <- 0
fixture$estimated_variance_2 <- 1; fixture$estimated_variance_5 <- 1
stopifnot(e$perceived_loss(fixture)$choice==2L,e$perceived_loss(fixture)$tied)
fixture$estimated_variance_5 <- .9
stopifnot(e$perceived_loss(fixture)$choice==5L)
# Failures remain explicit, rather than being mapped to an attractive loss.
bad <- obs; bad$D <- 0
stopifnot(!e$researcher_choice(bad)$decision$ok)
# Hand-calculated normalized gains, including negative/unclipped values.
expect_close(attainable_fraction(10,9,8),.5)
expect_close(attainable_fraction(10,12,8),-1)
expect_close(attainable_fraction(10,7,8),1.5)
stopifnot(is.na(attainable_fraction(10,10,10)))
small <- data.frame(oracle_loss_2=c(1,3,2),oracle_loss_5=c(2,1,4),oracle_choice=c(2,5,2))
perf <- researcher_performance(small,c(2L,2L,5L))
expect_close(perf$D_fixed,2); expect_close(perf$D_oracle,4/3)
expect_close(perf$D_feasible,8/3); expect_close(perf$fraction_captured,-1)
# Reproducibility covers choices, oracle evaluation, menu, and summaries.
pop <- population_estimands(p)
x <- lapply(1:3,function(r) researcher_replication(r,p,pop$theta_a))
y <- lapply(1:3,function(r) researcher_replication(r,p,pop$theta_a))
stopifnot(identical(x,y))
d <- do.call(rbind,lapply(x,`[[`,"pair"))
summary <- researcher_policy_summary(d)
expect_close(summary$total_bias,summary$systematic_distortion+summary$epsilon_distortion)
expect_close(summary$epsilon_distortion,summary$design_offset+summary$outcome_noise_distortion)
stopifnot(all(summary$fraction_captured[summary$attainable_gain>0]<=1+1e-12),
          all(is.na(summary$fraction_captured[summary$attainable_gain==0])))
# Exact fixed-Z Gaussian score covariance confirms the selection channel:
# Cov(eta5,gamma_hat)=0; Cov(eta2,gamma_hat)=-c5*Var(gamma_hat).
expect_close(d$evaluation_gamma_noise_cov_5,rep(0,nrow(d)))
expect_close(d$evaluation_gamma_noise_cov_2,-d$c5_GX*d$evaluation_oracle_gamma_variance)
# A separate R process reproduces without workspace state.
script <- tempfile(fileext=".R"); result <- tempfile(fileext=".rds")
writeLines(c('source("R/simulation_helpers.R"); load_project(); pop <- population_estimands()',
 'x <- lapply(1:3,function(r) researcher_replication(r,config$dgp,pop$theta_a))',
 sprintf('saveRDS(x,%s)',encodeString(result,quote='"'))),script)
status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)))
stopifnot(status==0L,identical(x,readRDS(result)))
unlink(c(script,result))
