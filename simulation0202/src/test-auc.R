set.seed(234)

calc_auc <- function(predictions, is_null){
    r <- rank(predictions)
    r1 <- sum(r[!is_null])
    n1 <- sum(!is_null)
    n2 <- sum(is_null)
    (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}

compare_calcs <- function(){
  n <- runif(1, 10, 1000)
  predictions <- rnorm(n)
  is_null <- rep(TRUE, n)
  is_null[seq(1, n / 2)] <- FALSE
  us <- calc_auc(predictions, is_null)
  pred <- ROCR::prediction(predictions, !is_null)
  them <- ROCR::performance(pred, "auc")@y.values[[1]]
  stopifnot(all.equal(us, them))
}

invisible(replicate(100, compare_calcs()))
