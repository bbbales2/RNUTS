#' Get multiple NUTS samples from a posterior
#'
#' @param num_samples
#' @param q0
#' @param h0
#' @param ham_system
#' @param integrator
#' @param max_treedepth
#' @param DEBUG
#'
#' @return
#' @export
sampleNuts <- function(numDraws, q0, h, ham, max_treedepth = 10) {
  draws = 10000
  qs = matrix(0, nrow = draws, ncol = length(q0))
  colnames(qs) = paste0("q", 1:length(q0))
  qs[1, ] = q0
  for(draw in 2:draws) {
    qs[draw, ] = oneSampleNuts(qs[draw - 1, ], h, ham, max_treedepth)$q
  }

  return(qs)
}

log_sum_exp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - log_sum_exp(x))
}

norm = function(v) {
  sqrt(sum(v^2))
}

min_uturn = function(q_plus, q_minus, p_plus, p_minus) {
  min(as.numeric(p_plus %*% (q_plus - q_minus)) / (norm(p_plus) * norm(q_plus - q_minus)),
      as.numeric(-p_minus %*% (q_minus - q_plus)) / (norm(p_minus) * norm(q_minus - q_plus)))
}

uturn = function(q_plus, q_minus, p_plus, p_minus) {
  no_uturn_forward <- as.numeric(p_plus %*% (q_plus - q_minus)) > 0
  no_uturn_backward <- as.numeric(-p_minus %*% (q_minus - q_plus)) > 0

  return(!(no_uturn_forward & no_uturn_backward))
}

#' @export
oneSampleNuts = function(q0, h, ham, max_treedepth = 10, DEBUG = FALSE, seed = NULL) {
  if(is.null(seed)) {
    seed = sample.int(.Machine$integer.max, 1)
  }

  set.seed(seed)

  p0 = ham$sampleMomentum()

  H0 = ham$H(q0, p0)

  # sample directions we'll go in ahead of time for easier debugging
  directions = c()
  depth_map = c(0)
  for(depth in 1:max_treedepth) {
    direction = sample(c(-1, 1), 1)
    directions = c(directions, direction)
    if(direction < 0) {
      depth_map = c(rep(depth, 2^max(0, depth - 1)), depth_map)
    } else {
      depth_map = c(depth_map, rep(depth, 2^max(0, depth - 1)))
    }
  }

  steps = list()
  qs = matrix(0, nrow = 2^max_treedepth, ncol = length(q0))
  ps = matrix(0, nrow = 2^max_treedepth, ncol = length(p0))
  log_pi = rep(0, 2^max_treedepth)
  i_first = which(depth_map == 0)[1]
  qs[i_first, ] = q0
  ps[i_first, ] = p0
  log_pi[i_first] = -H0
  i_left = i_first
  i_right = i_first
  uturn_detected = rep(0, max_treedepth)
  log_sum_pi_old = log_pi[i_first]
  i_old = i_first
  uturn_scores = c()
  uturn_left = c()
  uturn_right = c()
  for(depth in 1:max_treedepth) {
    if(directions[depth] < 0) {
      depth_steps = rev(which(depth_map == depth))
      i_left = tail(depth_steps, 1)
    } else {
      depth_steps = which(depth_map == depth)
      i_right = tail(depth_steps, 1)
    }
    steps[[depth]] = depth_steps

    for(i in depth_steps) {
      if(directions[depth] < 0) {
        i_prev = i + 1
      } else {
        i_prev = i - 1
      }

      z = leapfrog_step(qs[i_prev, ], ps[i_prev, ], h * directions[depth], ham)
      qs[i, ] = z$q
      ps[i, ] = z$p
      log_pi[i] = -ham$H(z$q, z$p)
    }

    uturn_detected_new_tree = FALSE
    tree_depth = log2(length(depth_steps))
    if(tree_depth > 0) {
      for(uturn_depth in tree_depth:1) {
        div_length = 2^(uturn_depth)
        starts = seq(1, length(depth_steps), div_length)
        ends = starts + div_length - 1

        if(directions[depth] < 0) {
          tmp = starts
          starts = ends
          ends = tmp
        }

        for(j in 1:length(starts)) {
          is_uturn = uturn(qs[depth_steps[starts[j]], ],
                             qs[depth_steps[ends[j]], ],
                             ps[depth_steps[starts[j]], ],
                             ps[depth_steps[ends[j]], ])
          uturn_scores = c(uturn_scores, min_uturn(qs[depth_steps[starts[j]], ],
                                                 qs[depth_steps[ends[j]], ],
                                                 ps[depth_steps[starts[j]], ],
                                                 ps[depth_steps[ends[j]], ]))
          uturn_left = c(uturn_left, depth_steps[ends[j]])
          uturn_right = c(uturn_right, depth_steps[starts[j]])

          if(is_uturn) {
            uturn_detected_new_tree = TRUE
            break
          }
        }

        if(uturn_detected_new_tree) {
          break
        }
      }
    }

    uturn_detected[depth] = uturn_detected_new_tree | uturn(qs[i_right, ], qs[i_left, ], ps[i_right, ], ps[i_left, ])
    uturn_scores = c(uturn_scores, min_uturn(qs[i_right, ], qs[i_left, ], ps[i_right, ], ps[i_left, ]))
    uturn_left = c(uturn_left, i_left)
    uturn_right = c(uturn_right, i_right)

    if(uturn_detected[depth]) {
      break
    }

    log_sum_pi_old = log_pi[i_first]
    log_sum_pi_new = log_sum_exp(log_pi[depth_steps])

    if(length(depth_steps) > 1) {
      i_new = sample(depth_steps, 1, prob = softmax(log_pi[depth_steps]))
    } else {
      i_new = depth_steps
    }

    i_old = sample(c(i_old, i_new), 1, prob = softmax(c(log_sum_pi_old, log_sum_pi_new)))
  }
  q = qs[i_old, ]

  if(DEBUG) {
    colnames(qs) = paste0("q", 1:length(q0))
    colnames(ps) = paste0("p", 1:length(p0))
    trajectory = bind_cols(c(tibble(log_pi = log_pi,
                       depth_map = depth_map),
                as_tibble(qs),
                as_tibble(ps)))[which(depth_map <= depth),]
    uturn_data = tibble(uturn_scores = uturn_scores,
                        uturn_left = uturn_left,
                        uturn_right = uturn_right) %>%
      mutate(depth = log2(uturn_right - uturn_left + 1))
    return(list(q = q,
                i = i_old,
                q0 = q0,
                h = h,
                max_treedepth = max_treedepth,
                seed = seed,
                uturn_data = uturn_data,
                trajectory = trajectory,
                directions = directions,
                uturn_detected = uturn_detected[1:depth],
                steps = steps[1:depth]))
  } else {
    return(list(q = q))
  }
}
