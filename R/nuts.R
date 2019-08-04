#' Generate multiple MCMC draws using multinomial nuts
#'
#' @param numDraws Number of draws to take
#' @inheritParams oneSampleNuts
#'
#' @return A matrix of draws. Columns are different parameters, rows are different draws
#' @export
sampleNuts <- function(numDraws, q0, h, ham, max_treedepth = 10) {
  qs = matrix(0, nrow = numDraws, ncol = length(q0))
  colnames(qs) = paste0("q", 1:length(q0))
  qs[1, ] = q0
  for(draw in 2:numDraws) {
    out = oneSampleNuts(qs[draw - 1, ], h, ham, max_treedepth)
    qs[draw, ] = out$q
  }

  return(qs)
}

log_sum_exp = function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax = function (x) {
  exp(x - log_sum_exp(x))
}

norm = function(v) {
  sqrt(sum(v^2))
}

uturn = function(q_plus, q_minus, p_plus, p_minus) {
  no_uturn_forward <- as.numeric(p_plus %*% (q_plus - q_minus)) > 0
  no_uturn_backward <- as.numeric(-p_minus %*% (q_minus - q_plus)) > 0

  return(!(no_uturn_forward & no_uturn_backward))
}

#' Generate a draw using Multinomial NUTS (https://arxiv.org/abs/1701.02434 with the original
#' Uturn criteria https://arxiv.org/abs/1111.4246)
#'
#' @param q0 current draw
#' @param h leapfrog stepsize
#' @param ham system hamiltonian (generated from hamwrapper library)
#' @param max_treedepth max treedepth
#' @param debug if true, then return a bunch of debugging info, otherwise return only the next draw
#' @param seed use the same seed to get the same behavior -- if NULL, make one up. Whatever is used is returned if debug == TRUE
#' @return if debug is false, list containing single element (q) containing next draw. If debug is true, a lot of stuff
#'
#' @export
oneSampleNuts = function(q0, h, ham, max_treedepth = 10, debug = FALSE, seed = NULL) {
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

        sorted_depth_steps = sort(depth_steps)

        for(j in 1:length(starts)) {
          is_uturn = uturn(qs[sorted_depth_steps[ends[j]], ],
                             qs[sorted_depth_steps[starts[j]], ],
                             ps[sorted_depth_steps[ends[j]], ],
                             ps[sorted_depth_steps[starts[j]], ])

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

    uturn_detected = uturn_detected_new_tree | uturn(qs[i_right, ], qs[i_left, ], ps[i_right, ], ps[i_left, ])

    if(uturn_detected) {
      break
    }

    log_sum_pi_new = log_sum_exp(log_pi[depth_steps])

    if(length(depth_steps) > 1) {
      i_new = sample(depth_steps, 1, prob = softmax(log_pi[depth_steps]))
    } else {
      i_new = depth_steps
    }

    p_new = min(1, exp(log_sum_pi_new - log_sum_pi_old))

    #print(paste("log_sum_pi_old", log_sum_pi_old))
    #print(paste("log_sum_pi_new", log_sum_pi_new))
    #print(paste("p_new", p_new))
    #print(softmax(c(log_sum_pi_old, log_sum_pi_new)))
    #print("--")

    i_old = sample(c(i_old, i_new), 1, prob = c(1 - p_new, p_new))

    log_sum_pi_old = log_sum_exp(c(log_sum_pi_old, log_sum_pi_new))
  }
  q = qs[i_old, ]

  if(debug) {
    colnames(qs) = paste0("q", 1:length(q0))
    colnames(ps) = paste0("p", 1:length(p0))
    trajectory = bind_cols(c(tibble(log_pi = log_pi,
                       depth_map = depth_map),
                as_tibble(qs),
                as_tibble(ps)))[which(depth_map <= depth),]
    return(list(q = q,
                i = i_old - min(which(depth_map <= depth)) + 1,
                q0 = q0,
                h = h,
                max_treedepth = max_treedepth,
                seed = seed,
                trajectory = trajectory,
                directions = directions,
                steps = steps))
  } else {
    return(list(q = q))
  }
}
