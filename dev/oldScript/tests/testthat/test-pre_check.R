test_that("test run_check", {
  expect_message(run_check(counts_input, group_list, tempdir()),"have done")
})
