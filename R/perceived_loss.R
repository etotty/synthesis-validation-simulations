# Feasible decision only. Every numeric input derives from the observed data.
perceived_loss <- function(information, precision="augmented", tie_tolerance=1e-10) {
  stopifnot(precision %in% c("augmented","HC3"), length(tie_tolerance)==1L,
            is.finite(tie_tolerance), tie_tolerance>=0)
  if (!information$ok) return(list(ok=FALSE, status=information$status, choice=NA_integer_))
  if (precision=="augmented") {
    v2 <- information$estimated_variance_2; v5 <- information$estimated_variance_5
  } else { v2 <- information$hc3_variance_2; v5 <- information$hc3_variance_5 }
  penalty <- information$estimated_bias_5^2
  if (any(!is.finite(c(v2,v5,penalty))) || min(v2,v5,penalty)<0)
    return(list(ok=FALSE, status="invalid_precision_or_penalty", choice=NA_integer_))
  l2 <- v2; l5 <- v5+penalty
  # Within an absolute tolerance, prefer the flexible reference analysis 2.
  list(ok=TRUE, status="ok", precision=precision, perceived_loss_2=l2,
       perceived_loss_5=l5, perceived_bias_penalty_2=0,
       perceived_bias_penalty_5=penalty, variance_proxy_2=v2, variance_proxy_5=v5,
       choice=if (l2 <= l5+tie_tolerance) 2L else 5L,
       tied=abs(l2-l5)<=tie_tolerance)
}

researcher_choice <- function(data, precision="augmented", tie_tolerance=1e-10) {
  info <- researcher_information(data)
  list(information=info, decision=perceived_loss(info, precision, tie_tolerance))
}
