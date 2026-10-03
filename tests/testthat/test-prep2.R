library(testthat)


test_that("prep2 reproduces every decision-tree leaf", {
  expect_equal(prep2(8, 65)$category, "Excellent")
  expect_equal(prep2(6, 82)$category, "Good")
  expect_equal(prep2(3, 60, "present")$category, "Good")
  expect_equal(prep2(2, 70, "absent", 5)$category, "Limited")
  expect_equal(prep2(1, 74, "absent", 12)$category, "Poor")
})


test_that("prep2 locks cutoff boundary semantics", {
  at_safe <- prep2(5, 80)
  at_nihss <- prep2(4, 60, mep_status = FALSE, nihss = 7)

  expect_equal(at_safe$category, "Good")
  expect_equal(at_safe$pathway, c("SAFE>=5", "age>=80"))
  expect_equal(at_nihss$category, "Poor")
  expect_equal(
    at_nihss$pathway,
    c("SAFE<5", "MEP-", "NIHSS>=7")
  )
})


test_that("prep2 accepts documented MEP aliases", {
  expect_equal(prep2(3, 60, TRUE)$category, "Good")
  expect_equal(prep2(3, 60, "MEP+")$category, "Good")
  expect_equal(prep2(3, 60, "PRESENT")$category, "Good")
  expect_equal(prep2(3, 60, FALSE, 4)$category, "Limited")
  expect_equal(prep2(3, 60, "MEP-", 10)$category, "Poor")
})


test_that("prep2 requires branch-specific inputs", {
  expect_error(prep2(2, 60), "MEP status")
  expect_error(prep2(2, 60, "absent"), "NIHSS")
  expect_silent(prep2(8, 60))
})


test_that("prep2 validates clinical ranges", {
  expect_error(prep2(11, 60), "safe_score")
  expect_error(prep2(2.5, 60), "safe_score")
  expect_error(prep2(5, 0), "age")
  expect_error(prep2(2, 60, "unknown"), "mep_status")
  expect_error(prep2(2, 60, "absent", 43), "nihss")
})


test_that("prep2 result prints and echoes normalized inputs", {
  result <- prep2(2, 70, "MEP-", 5)

  expect_s3_class(result, "prep2_result")
  expect_equal(result$inputs$mep_status, "absent")
  expect_match(result$description, "Modest")
  expect_output(print(result), "PREP2 prognosis: Limited")
})
