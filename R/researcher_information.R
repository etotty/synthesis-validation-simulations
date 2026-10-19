# Observable-only module. No configuration, DGP parameters, propensity scores,
# conditional means or population estimands are accepted or looked up here.
# Source with perceived_loss.R in an environment containing only base/stats.
observed_ols <- function(W, y, contrast) {
  if (nrow(W) <= ncol(W)) return(list(ok=FALSE, status="insufficient_rows"))
  qrW <- qr(W, tol=1e-10)
  if (qrW$rank < ncol(W)) return(list(ok=FALSE, status="rank_deficient"))
  Q <- qr.Q(qrW)
  h <- drop(Q %*% backsolve(qr.R(qrW), contrast[qrW$pivot], transpose=TRUE))
  residual <- qr.resid(qrW, y)
  leverage <- rowSums(Q^2)
  hc3 <- if (any(1-leverage <= 1e-12)) NA_real_ else sum((h*residual/(1-leverage))^2)
  list(ok=TRUE, status="ok", coefficients=qr.coef(qrW,y), h=h,
       estimate=sum(h*y), variance_hc3=hc3, weight_norm=sum(h^2),
       residual_variance=sum(residual^2)/(nrow(W)-ncol(W)),
       residual=residual, leverage=leverage, W=W, qr=qrW)
}

researcher_information <- function(data) {
  allowed <- c("Y", "X", "D", "G")
  if (!is.data.frame(data) || !setequal(names(data), allowed) || ncol(data)!=4L)
    stop("Researcher input must contain exactly Y, X, D, G; project evaluation data first")
  stopifnot(all(vapply(data, is.numeric, logical(1))), all(is.finite(as.matrix(data))),
            all(data$G %in% 0:1), all(data$D %in% 0:1))
  n <- nrow(data)
  W5 <- cbind(intercept=rep(1,n), D=data$D, X=data$X, X2=data$X^2,
              G=data$G, DG=data$D*data$G)
  c5 <- c(0,1,0,0,0,1)
  sub <- which(data$G==1)
  W2 <- W5[sub, c("intercept","D","X","X2"), drop=FALSE]
  Waug <- cbind(W5, GX=data$G*data$X)
  fit2 <- observed_ols(W2, data$Y[sub], c(0,1,0,0))
  fit5 <- observed_ols(W5, data$Y, c5)
  # Augmented contrast is the subgroup treatment effect; GX is a diagnostic,
  # not a sixth selectable analysis.
  aug <- observed_ols(Waug, data$Y, c(c5,0))
  fits <- list(subgroup=fit2, pooled=fit5, augmented=aug)
  failed <- !vapply(fits, function(x) x$ok, logical(1))
  if (any(failed)) return(list(ok=FALSE,
    status=paste(names(fits)[failed], vapply(fits[failed], function(x) x$status, character(1)), collapse=";")))
  gamma_contrast <- c(rep(0,6),1)
  gamma_h <- drop(qr.Q(aug$qr) %*% backsolve(qr.R(aug$qr),
                         gamma_contrast[aug$qr$pivot], transpose=TRUE))
  h2 <- numeric(n); h2[sub] <- fit2$h
  sigma2_hat <- aug$residual_variance
  gamma_hat <- unname(aug$coefficients["GX"])
  c5_GX <- sum(fit5$h * data$G * data$X)
  list(ok=TRUE, status="ok", gamma_hat=gamma_hat,
       gamma_se=sqrt(sigma2_hat*sum(gamma_h^2)),
       sigma2_hat=sigma2_hat, c5_GX=c5_GX, estimated_bias_5=gamma_hat*c5_GX,
       estimated_variance_2=sigma2_hat*fit2$weight_norm,
       estimated_variance_5=sigma2_hat*fit5$weight_norm,
       hc3_variance_2=fit2$variance_hc3, hc3_variance_5=fit5$variance_hc3,
       T_2=fit2$estimate, T_5=fit5$estimate, T_augmented=aug$estimate,
       h2=h2, h5=fit5$h, gamma_h=gamma_h,
       n=n, n_subgroup=length(sub), n_subgroup_treated=sum(data$D[sub]),
       n_subgroup_control=sum(data$D[sub]==0),
       max_leverage_2=max(fit2$leverage), max_leverage_5=max(fit5$leverage))
}
