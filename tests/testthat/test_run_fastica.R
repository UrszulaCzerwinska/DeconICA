library(fastICA)

test_that("fastICA output dimensions", {
  S <- matrix(runif(10000), 5000, 2)
  A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
  X <- data.frame(S %*% A)
  test <-
    deconica::run_fastica(
      X,
      row.center = TRUE,
      n.comp = 2,
      optimal = FALSE,
      with.names = FALSE
    )
  expect_equal(ncol(test$S), 2)
  expect_equal(nrow(test$S), nrow(X))
  expect_equal(nrow(test$A), 2)
  expect_equal(ncol(test$A), ncol(X))
})

test_that("fastICA testing optimal for limited number of samples", {
  set.seed(10)
  S <- matrix(runif(10000), 5000, 70)
  A <- matrix(c(1, 1, -1, 3), 70, 90, byrow = TRUE)
  X <- data.frame(S %*% A)
  test <-
    deconica::run_fastica(
      X,
      row.center = TRUE,
      n.comp = 10,
      optimal = TRUE,
      with.names = FALSE
    )
  expect_equal(ncol(test$S), 2)
  expect_equal(nrow(test$S), nrow(X))
})

test_that("fastICA testing optimal for big number of samples", {
  set.seed(10)
  S <- matrix(runif(10000), 5000, 70)
  A <- matrix(c(1, 1, -1, 3), 70, 200, byrow = TRUE)
  X <- data.frame(S %*% A)
  test <-
    deconica::run_fastica(
      X,
      row.center = TRUE,
      n.comp = 10,
      optimal = TRUE,
      with.names = FALSE
    )
  expect_equal(ncol(test$S), 100)
  expect_equal(nrow(test$S), nrow(X))
})

test_that("fastICA testing names", {
  set.seed(10)
  S <- matrix(runif(10000), 5000, 70)
  A <- matrix(c(1, 1, -1, 3), 70, 200, byrow = TRUE)
  X <- data.frame(S %*% A)
  names <- paste("A", 1:nrow(X), sep = "")
  X <- cbind(names, X)
  test <-
    deconica::run_fastica(
      X,
      row.center = TRUE,
      n.comp = 10,
      optimal = TRUE,
      with.names = TRUE
    )
  expect_equal(test$names, names)
  expect_equal(nrow(test$S), nrow(X))
})

test_that("removing duplicates", {
  set.seed(10)
  X <- matrix(runif(3500), 50, 70)
  names <- paste("A", sample(1:10, nrow(X), replace = TRUE), sep = "")
  X <- data.frame(names, X)
  X.2 <- .remove_duplicates(X)
  expect_equal(nrow(X.2), length(unique(names)))
})
