library(rstan)
library(hamwrapper)

test_that("we can generate a sample", {
  testModelFile = system.file("extdata", "normal.stan", package = "RNUTS")

  ham = createHamiltonianSystem(testModelFile, list(N = 5))

  q0 = rnorm(5)
  out = oneSampleNuts(q0, 0.1, ham)
  expect_equal(length(out$q), 5)
})

test_that("we're roughly sampling a normal correctly", {
  testModelFile = system.file("extdata", "normal.stan", package = "RNUTS")

  ham = createHamiltonianSystem(testModelFile, list(N = 1))

  q0 = rnorm(1)
  N = 10000
  qs = sampleNuts(N, q0, 1.0, ham)

  Neff = rstan::monitor(array(qs[,1], dim = c(nrow(qs), 1, 1)), warmup = 0, print = FALSE)$n_eff
  expect_gt(Neff, 1000)

  K = round(2 * Neff^0.4)

  quantiles = qnorm((0:K) / K)
  counts = hist(qs[,1], quantiles, plot = FALSE)$counts * Neff / N
  expected = rep(Neff / K, K)
  score = sum((counts - expected)^2 / expected)

  expect_lt(score, qchisq(1e-6, df = K - 1, lower.tail = FALSE))
})
