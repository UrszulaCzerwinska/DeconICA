test_that("correlate_metagenes dimensions output", {
  res_run_ica <- deconica::run_fastica (deconica::Example_ds, optimal = FALSE, n.comp = 5,
             with.names = TRUE)
   test <-
    deconica::correlate_metagenes(
      S = res_run_ica$S,
      gene.names = res_run_ica$names)
  expect_equal(nrow(test$n),ncol(res_run_ica$S) )
  expect_equal(ncol(test$n), length(data.list))
  expect_equal(length(test), 3)
})

test_that("correlate_metagenes strong threshold", {
  res_run_ica <- deconica::run_fastica (deconica::Example_ds, optimal = FALSE, n.comp = 5,
                                        with.names = TRUE)
  test <-
    deconica::correlate_metagenes(
      S = res_run_ica$S,
      gene.names = res_run_ica$names, threshold = 4)
  expect_equal(nrow(test$n),0 )
  expect_equal(ncol(test$n),11)
  expect_equal(length(test$n), 0)
  expect_equal(length(test), 3)
})

test_that("output assign_metagenes", {
  set.seed(299)
  res_run_ica <- deconica::run_fastica (deconica::Example_ds, optimal = FALSE, n.comp = 5,
                                        with.names = TRUE)
  res_corr_matagenes <-
    deconica::correlate_metagenes(
      S = res_run_ica$S,
      gene.names = res_run_ica$names)
  test <- deconica::assign_metagenes(res_corr_matagenes$r)
  expect_equal(ncol(test), 2L)
  expect_true(is.factor(test[[2]]))
  expect_true(is.factor(test[[1]]))
  expect_equal(names(test), c("metagene","IC"))
})

test_that("output identify_immune_ic", {
  set.seed(299)
  res_run_ica <- deconica::run_fastica (deconica::Example_ds, optimal = FALSE, n.comp = 5,
                                        with.names = TRUE)
  res_corr_matagenes <-
    deconica::correlate_metagenes(
      S = res_run_ica$S,
      gene.names = res_run_ica$names)
  res_assign_metagenes <- deconica::assign_metagenes(res_corr_matagenes$r)
  test <- deconica::identify_immune_ic(res_corr_matagenes$r[,"M8_IMMUNE"], res_assign_metagenes[[2]])
  test2 <- deconica::identify_immune_ic(res_corr_matagenes$r[,"M8_IMMUNE"], "")
  expect_equal(length(test), 0)
  expect_equal(length(test2), 2)
  expect_equal(typeof(test2), typeof(test))
})


