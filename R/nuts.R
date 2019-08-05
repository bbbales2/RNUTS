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

# This is the uturn criteria defined in Eq. 3 of A.4.2 https://arxiv.org/pdf/1701.02434.pdf or
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

  # directions is a length max_treedepth vector that maps each treedepth
  #  to an integration direction (left or right)
  directions = c()
  # depth_map will be a vector of length 2^max_treedepth that maps each of
  #  the possibly 2^max_treedepth points to a treedepth
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

  # Steps is a list that maps treedepth to which leapfrog steps were
  #  computed in that treedepth (kinda the opposite of depth_map)
  steps = list()
  # qs stores our positions
  qs = matrix(0, nrow = 2^max_treedepth, ncol = length(q0))
  # ps stores our momentums
  ps = matrix(0, nrow = 2^max_treedepth, ncol = length(p0))
  # log_pi defined in section A.2.3 of https://arxiv.org/abs/1701.02434
  log_pi = rep(0, 2^max_treedepth)
  # index of initial state
  i_first = which(depth_map == 0)[1]
  qs[i_first, ] = q0
  ps[i_first, ] = p0
  log_pi[i_first] = -H0
  # i_left and i_right are indices that track the leftmost and rightmost
  #  states of the integrated trajectory
  i_left = i_first
  i_right = i_first
  # log_sum_pi_old is the log of the sum of the pis (of log_pi) for the
  #  tree processed so far
  log_sum_pi_old = log_pi[i_first]
  # i_old will be the sample chosen (the sample from T(z' | told) in section
  #  A.3.1 of https://arxiv.org/abs/1701.02434)
  i_old = i_first
  # For trees of increasing treedepth
  for(depth in 1:max_treedepth) {
    # Figure out what leapfrog steps we need to compute. If integrating in the
    #  positive direction update the index that points at the right side of the
    #  trajectory. If integrating in the negative direction update the index pointing
    #  to the left side of the trajectory.
    if(directions[depth] < 0) {
      depth_steps = rev(which(depth_map == depth))
      i_left = tail(depth_steps, 1)
    } else {
      depth_steps = which(depth_map == depth)
      i_right = tail(depth_steps, 1)
    }
    steps[[depth]] = depth_steps

    # Actually do the integrationg
    for(i in depth_steps) {
      if(directions[depth] < 0) {
        i_prev = i + 1
      } else {
        i_prev = i - 1
      }

      # This doesn't take advantage of the leapfroggy-nature of leapfrog.
      z = leapfrog_step(qs[i_prev, ], ps[i_prev, ], h * directions[depth], ham)
      qs[i, ] = z$q
      ps[i, ] = z$p
      log_pi[i] = -ham$H(z$q, z$p)
    }

    # What we're doing here is generating a trajectory composed of a number of leapfrog states.
    # We apply a set of comparisons on this trajectory that can be organized to look like a binary tree.
    # Sometimes I say trajectory instead of tree. Each piece of the tree in the comparison corresponds
    #  to a subset of the trajectory. When I say trajectory I'm referring to a piece of the trajectory that
    #  also corresponds to some sorta subtree in the comparisons.
    #
    # This is probably confusing but what I want to communicate is trajectory and tree are very related
    #  but maybe technically not the same thing.
    #
    # Detect U-turns in newly integrated subtree
    uturn_detected_new_tree = FALSE
    tree_depth = log2(length(depth_steps))
    if(tree_depth > 0) {
      # Start at root of new subtree and work down to leaves
      for(uturn_depth in tree_depth:1) {
        # The root of the comparison tree compares the leftmost to the rightmost states of the new
        #  part of the trajectory.
        #  The next level down in the tree of comparisons cuts that trajectory in two and compares
        #  the leftmost and rightmost elements of those smaller trajectories.
        #  Etc. Etc.
        div_length = 2^(uturn_depth)

        # Starts are relative indices pointing to the leftmost state of each comparison to be done
        starts = seq(1, length(depth_steps), div_length)
        # Ends are relative indices pointing to the rightmost state for each comparison to be done
        ends = starts + div_length - 1

        # Starts and ends are relative because they point to the ith leapfrog step in a sub-trajectory
        #  of size 2^tree_depth which needs to be mapped to the global index of qs
        #  (which is defined in the range 1:2^max_treedepth)
        #
        # The sort is necessary because depth_steps is sorted in order of leapfrog steps taken
        #  which might be backwards in time (so decreasing instead of increasing)
        #
        # This sort is important because we need to keep track of what is left and what is right
        #  in the trajectory so that we do the right comparisons
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

    # Merging the two trees requires one more uturn check from the overall left to right states
    uturn_detected = uturn_detected_new_tree | uturn(qs[i_right, ], qs[i_left, ], ps[i_right, ], ps[i_left, ])

    if(uturn_detected) {
      break
    }

    # log of the sum of pi (A.3.1 of https://arxiv.org/abs/1701.02434) of the new subtree
    log_sum_pi_new = log_sum_exp(log_pi[depth_steps])

    # sample from the new subtree according to the equation in A.2.1 in https://arxiv.org/abs/1701.02434
    #  (near the end of that section)
    if(length(depth_steps) > 1) {
      i_new = sample(depth_steps, 1, prob = softmax(log_pi[depth_steps]))
    } else {
      i_new = depth_steps
    }

    # Pick between the samples generated from the new and old subtrees using the biased progressive sampling in
    #  A.3.2 of https://arxiv.org/abs/1701.02434
    p_new = min(1, exp(log_sum_pi_new - log_sum_pi_old))
    i_old = sample(c(i_old, i_new), 1, prob = c(1 - p_new, p_new))

    # Update log of sum of pi of overall tree
    log_sum_pi_old = log_sum_exp(c(log_sum_pi_old, log_sum_pi_new))
  }
  # Get the final sample
  q = qs[i_old, ]

  if(debug) {
    colnames(qs) = paste0("q", 1:length(q0))
    colnames(ps) = paste0("p", 1:length(p0))
    # For a system with N parameters, this tibble will have
    #  2N + 2 columns. The first column (log_pi) stores log of pi (A.2.3 of https://arxiv.org/abs/1701.02434)
    #  The second column (depth_map) stores at what treedepth that state was added to the trajectory
    #  The next N columns are the positions (all starting with q)
    #  The next N columns are the momentums (all starting with p)
    trajectory = bind_cols(c(tibble(log_pi = log_pi,
                       depth_map = depth_map),
                as_tibble(qs),
                as_tibble(ps)))[which(depth_map <= depth),]
    return(list(q = q, # new sample
                i = i_old - min(which(depth_map <= depth)) + 1, # q is the ith row of trajectory
                q0 = q0,
                h = h,
                max_treedepth = max_treedepth,
                seed = seed, # seed used by random number generator, for reproducibility
                trajectory = trajectory, # tibble containing details of trajectory
                directions = directions)) # direction integrated in each subtree
  } else {
    return(list(q = q))
  }
}
