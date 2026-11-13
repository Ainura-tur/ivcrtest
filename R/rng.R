# Random-number state helpers.
#
# Several functions take a `seed` argument so that their results are
# reproducible. A bare set.seed() inside a function overwrites .Random.seed in
# the user's global environment, which CRAN policy forbids and which makes
# later draws in the session repeat. .cr_local_seed() instead sets the seed for
# the lifetime of the calling function and restores the caller's state when
# that function exits, normally or through an error.

# Set `seed` until the calling function exits. With seed = NULL nothing is set
# or restored, so the session's stream is used and advanced as usual.
.cr_local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(NULL))
  old <- .cr_rng_state()
  restore <- as.call(list(function() .cr_restore_rng_state(old)))
  do.call(base::on.exit, list(restore, TRUE, FALSE), envir = envir)
  set.seed(seed)
  invisible(old)
}

# The current state, or NULL when the session has not used the RNG yet.
.cr_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.cr_restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}
