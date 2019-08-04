library(tidyverse)
library(rstan)
library(ggplot2)
library(hamwrapper)

# Use the test model included in the package
testModelFile = system.file("extdata", "normal.stan", package = "RNUTS")

# Load up a model and generate a bunch of samples
#   N here is the dimension of the std_normal we'll be sampling
ham = createHamiltonianSystem(testModelFile, list(N = 1))
q0 = rnorm(1)
N = 1000
qs = sampleNuts(N, q0, 0.5, ham)
rstan::monitor(array(qs[,1], dim = c(nrow(qs), 1, 1)), warmup = 0, print = FALSE)$n_eff

# Load up a different model and build a NUTS trajectory with debug info turned on
#   Scroll down to the end of the file to see a plot diagnosing the U-turn
ham = createHamiltonianSystem(testModelFile, list(N = 2))
q0 = rnorm(2)
out = oneSampleNuts(q0, 0.1, ham, debug = TRUE, seed = 633273650)

# Compute minimum of the two uturn criteria
min_uturn = function(q_plus, q_minus, p_plus, p_minus) {
  min(as.numeric(p_plus %*% (q_plus - q_minus)) / (norm(p_plus) * norm(q_plus - q_minus)),
      as.numeric(-p_minus %*% (q_minus - q_plus)) / (norm(p_minus) * norm(q_minus - q_plus)))
}

# This checks all the uturn criteria after the fact. I took this out of the
# sampler implementation so that the implementation would only have code that
# is necessary for
#
#uturn_detected = FALSE
tree_depth = log2(nrow(out$trajectory))
qs = out$trajectory %>% select(starts_with("q")) %>% as.matrix
ps = out$trajectory %>% select(starts_with("p")) %>% as.matrix
is_uturn = c()
uturn_scores = c()
uturn_left = c()
uturn_right = c()
if(tree_depth > 0) {
  for(uturn_depth in tree_depth:1) {
    div_length = 2^(uturn_depth)
    starts = seq(1, nrow(out$trajectory), div_length)
    ends = starts + div_length - 1

    for(j in 1:length(starts)) {
      is_uturn = c(is_uturn, uturn(qs[ends[j], ],
                                   qs[starts[j], ],
                                   ps[ends[j], ],
                                   ps[starts[j], ]))

      uturn_scores = c(uturn_scores, min_uturn(qs[ends[j], ],
                                               qs[starts[j], ],
                                               ps[ends[j], ],
                                               ps[starts[j], ]))
      uturn_left = c(uturn_left, starts[j])
      uturn_right = c(uturn_right, ends[j])

      # Checking all uturns for now
      #if(is_uturn) {
      #  uturn_detected = TRUE
      #  break
      #}
    }

    #if(uturn_detected) {
    #  break
    #}
  }
}
uturn = tibble(scores = uturn_scores,
               left = uturn_left,
               right = uturn_right)

uturn %>%
  ggplot(aes(right - left, scores)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_log10() +
  xlab("Trajectory length") +
  ggtitle("score below the red line (zero) indicates uturn")
