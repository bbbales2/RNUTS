#' @export
leapfrog_step <- function(q0, p0, h, ham) {
  GradU = function(q) ham$gradU(q)
  M = ham$M

  # leapfrog step
  p_half = p0 - (h / 2) * GradU(q0)
  q1 = q0 + h * solve(M, p_half)
  p1 = p_half - (h / 2) * GradU(q1)

  # return new state
  return(list(q = q1, p = p1))
}
