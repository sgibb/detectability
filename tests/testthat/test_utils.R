context("utils")

test_that(".msg", {
  expect_message(detectability:::.msg(TRUE, "foobar"), "foobar")
  expect_message(detectability:::.msg(TRUE, "foo", "bar"), "foobar")
  expect_silent(detectability:::.msg(FALSE, "foobar"))
})
