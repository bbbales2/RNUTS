library(tidyverse)
library(rstan)
library(ggplot2)
library(hamwrapper)
library(RNUTS)

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
out = oneSampleNuts(q0, 0.1, ham, debug = TRUE)

# Compute minimum of the two uturn criteria
min_uturn = function(q_plus, q_minus, p_plus, p_minus) {
  min(as.numeric(p_plus %*% (q_plus - q_minus)) / (norm(p_plus) * norm(q_plus - q_minus)),
      as.numeric(-p_minus %*% (q_minus - q_plus)) / (norm(p_minus) * norm(q_minus - q_plus)))
}

# This checks all the uturn criteria after the fact. I took this out of the
# sampler implementation so that the implementation so that would only have code
# absolutely necessary to sampling. Not sure if it's better out or in, honestly.
#
#uturn_detected = FALSE
tree_depth = log2(nrow(out$trajectory))
qs = out$trajectory %>% select(starts_with("q")) %>% as.matrix
ps = out$trajectory %>% select(starts_with("p")) %>% as.matrix
is_uturn = c()
uturn_scores = c()
uturn_left = c()
uturn_right = c()
uturn_depths = c()
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
      uturn_depths = c(uturn_depths, uturn_depth)

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

# The first uturn we hit will be the uturn that happens earliest in the last
#  trajectory expansion. To figure out which uturn this is, we need to figure
#  out which point is the first step in the last trajectory expansion
depth_map = out$trajectory %>% pull(depth_map)
last_integration_steps = which(depth_map == max(depth_map))
if(min(last_integration_steps) == 1) {
  initial_last_integration_step = max(last_integration_steps)
} else {
  initial_last_integration_step = min(last_integration_steps)
}

# Mark uturn checks that touch the invalid tree as invalid
uturn = tibble(scores = uturn_scores,
               left = uturn_left,
               right = uturn_right) %>%
  mutate(valid = out$trajectory[left,] %>% pull(valid) & out$trajectory[right,] %>% pull(valid))

# Figure out the first uturn (the first one we hit during the last expansion)
#  that would have made the tree invalid
invalid_uturn = uturn %>%
  filter(valid == FALSE,
         scores < 0.0) %>%
  mutate(distance = pmax(abs(left - initial_last_integration_step),
                         abs(right - initial_last_integration_step))) %>%
  top_n(1, -distance) %>%
  select(-distance)

# Plot the uturn scores as a function of distance
uturn %>%
  filter(valid == TRUE) %>%
  bind_rows(invalid_uturn) %>%
  ggplot(aes(right - left, scores)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_log10() +
  xlab("Trajectory length") +
  ylab("Uturn score\nIt's just the minimum of the two checks done") +
  ggtitle("Score below the red line (zero) indicates uturn")

# This is some stuff to get the color scale set up
min_scores = min(pull(uturn, scores))
max_scores = max(pull(uturn, scores))
rescale_to_0_1 = function(v) {
  (v - min(v)) / (max(v) - min(v))
}

# Plot all the uturn checks in tree form
uturn %>%
  mutate(depth = right - left + 1) %>%
  ggplot(aes()) +
  geom_curve(aes(x = left, xend = (left + right) / 2.0, y = 0, yend = depth, color = scores, linetype = valid), curvature = -0.25) +
  geom_curve(aes(x = (left + right) / 2.0, xend = right, y = depth, yend = 0, color = scores, linetype = valid), curvature = -0.25) +
  scale_color_gradientn(values = rescale_to_0_1(c(min_scores, -1e-8, 1e-8, max_scores)), colors = c("red", "red", "green", "blue")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  xlab("Point in trajectory") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle("Uturn checks\nA score less than zero is a uturn.\nMake sure to look at the discontinuity at 0.0 in the color scale")
