test_that("An example for a simulation study", {
  result=GDILM_SEIRS_Sim_Par_Est(3,3,8,30,0.7, 0.5, -1, 2.5, 0,30, 50,0.5,0.5, 2, 3, 10, 2)
  expect_type(result, "list")
})
