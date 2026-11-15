#' Re-score a population with replication, and pick the best members by their mean
#'
#' DEoptim evaluates each surviving member only once, so on a noisy objective the values it holds are
#' partly luck: a lucky low draw survives and is reported as "best". Re-scoring the final population
#' several times and ranking by the mean sees through that. In eight FireSense fits (2026-09-17),
#' DEoptim's best member ranked 1st to 6th of 60 by its replicated mean, and the ledger's "5 best" were
#' five copies of that one member.
#'
#' `rescorePopulation()` evaluates every row of `pop` `reps` times; each evaluation has its own seed
#' (`seed` + its index), so the result is reproducible. With a cluster the evaluations run on its
#' workers, which must already hold whatever `fn` reads from their global environment (as for DEoptim).
#'
#' @param pop A matrix, one parameter set per row, with the parameter names as column names.
#' @param fn The objective function; called as `fn(par, ...)` with `fnArgs`.
#' @param reps Evaluations per member.
#' @param cl A cluster from [parallel::makeCluster()], or `NULL` to evaluate here.
#' @param seed Base seed.
#' @param fnArgs Named list of further arguments to `fn`.
#'
#' @return `rescorePopulation()`: a `data.table` with `member` (row of `pop`), `rep` and `value`.
#' @export
#' @rdname replicatedBest
rescorePopulation <- function(pop, fn, reps = 10L, cl = NULL, seed = 1L, fnArgs = list()) {
  pop <- as.matrix(pop)
  jobs <- data.table::CJ(member = seq_len(NROW(pop)), rep = seq_len(reps))
  jobList <- lapply(seq_len(NROW(jobs)), function(i)
    list(par = pop[jobs$member[i], ], seed = as.integer(seed) + i))
  vals <- .runJobs(jobList, fn, fnArgs, cl)
  jobs[, value := as.numeric(unlist(vals))]
  jobs[]
}

#' Every member's objective terms, on common random numbers
#'
#' Evaluates every row of `pop` `reps` times with `fn(returnTerms = TRUE)`: replicate `r` of EVERY member uses
#' seed `seed + r`, so members are compared on the same random numbers ([rescorePopulation()] gives each
#' evaluation its own seed and keeps only the total). The terms can then be re-weighted afterwards without
#' simulating again.
#'
#' @inheritParams rescorePopulation
#' @param reps Evaluations of each member, `rep` 1 to `reps`.
#'
#' @return a `data.table`, one row per member x rep: `member` (row of `pop`), `rep`, and what `fn` returns with
#'   `returnTerms = TRUE`. For [.objfunSpreadFit()]: `objective`; `adTest`, `pitAD`, `yearArea` and `areaDist`
#'   as the objective added them (WEIGHTED); `crpsFire`, `crpsTotalArea`, `crpsTop10`, `crpsAnnualArea` and
#'   `crpsTopShare` unweighted; `nFires`; then, to re-weight, `raw_adTest`, `raw_pitAD`, `raw_yearArea` and
#'   `raw_areaDist` (before their weights) and the weight of every term, `w_<term>`. `objective` is the sum over
#'   terms of weight x unweighted value. A term the objective did not use is `NA` (`adTest`, `pitAD`, `yearArea` and `areaDist` when their tests are off); with `returnTerms`, `crpsAnnualArea` and `crpsTopShare` are computed even at weight 0, so they can be weighted up afterwards.
#' @export
scorePopulationTerms <- function(pop, fn = .objfunSpreadFit, reps = 1L, cl = NULL, seed = 1L, fnArgs = list()) {
  pop <- as.matrix(pop)
  fnArgs$returnTerms <- TRUE
  jobs <- data.table::CJ(member = seq_len(NROW(pop)), rep = seq_len(reps))
  jobList <- lapply(seq_len(NROW(jobs)), function(i)
    list(par = pop[jobs$member[i], ], seed = as.integer(seed) + jobs$rep[i]))
  vals <- .runJobs(jobList, fn, fnArgs, cl)
  cbind(jobs, data.table::as.data.table(do.call(rbind, vals)))
}

## Runs each job of `jobList` (`par`, `seed`), here or on the workers of `cl`
.runJobs <- function(jobList, fn, fnArgs, cl = NULL) {
  if (is.null(cl)) {
    lapply(jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  } else {
    parallel::clusterApplyLB(cl, jobList, .rescoreJob, fn = fn, fnArgs = fnArgs)
  }
}

## One evaluation, at namespace level so the workers receive it by reference, not with a caller's frame
.rescoreJob <- function(job, fn, fnArgs) {
  set.seed(job$seed)
  do.call(fn, c(list(par = job$par), fnArgs))
}

#' @param scores The output of `rescorePopulation()` for `pop`.
#' @param n How many members to return.
#'
#' @return `bestByReplicatedMean()`: a list with `params` (a `data.table`, one row per member, best
#'   first), `objFunVal` (their replicated means), `objFunValSD` and `member` (their rows in `pop`).
#'   Identical rows of `pop` count once, so the members returned are distinct.
#' @export
#' @rdname replicatedBest
bestByReplicatedMean <- function(pop, scores, n = 5L) {
  pop <- as.matrix(pop)
  keys <- apply(pop, 1, function(r) paste(sprintf("%a", r), collapse = ","))
  sc <- data.table::as.data.table(scores)
  ## identical rows of `pop` are one member: their evaluations are pooled
  sc[, key := keys[member]]
  pooled <- sc[, .(mean = mean(value), sd = stats::sd(value), member = min(member)), by = key]
  data.table::setorderv(pooled, "mean")
  top <- utils::head(pooled, n)
  list(params = data.table::as.data.table(pop[top$member, , drop = FALSE]),
       objFunVal = top$mean, objFunValSD = top$sd, member = top$member)
}

utils::globalVariables(c("member", "value", "key"))
