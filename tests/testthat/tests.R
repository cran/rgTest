testthat::test_that("test of wei_test", {
  d=200
  vmu2 = rep(1.7/sqrt(d),d) 
  num1 = 100
  num2 = 100
  y1a = matrix(0,num1,d)
  y2a = matrix(0,num2,d)
  
  for (i in 1:num1) {
    y1a[i,] = rlnorm(d, meanlog = 0, sdlog = 1)
  }
  for (i in 1:(num2)) {
    y2a[i,] = rlnorm(d, meanlog = vmu2, sdlog = 1)
  }
  y = rbind(y1a, y2a)
  dist = getdis(y)
  #no data, no distance matrix, no edge matrix
  testthat::expect_error(rg.test(n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 30))
  #no number of observations in the sample
  testthat::expect_error(rg.test(dis = dist, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 30))
  #wrong test.type
  testthat::expect_error(rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 30, test.type = list("so")))
  #wrong number of permutations
  testthat::expect_error(rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 'a'))
  #wrong number of permutations
  testthat::expect_message(rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 10001))
  #expect no error
  testthat::expect_no_error(rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 0))
  testthat::expect_no_error(rg.test(dis = dist, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 50))
  testthat::expect_no_error(rg.test(data.X = y1a, data.Y = y2a, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 50))
  E = kmst(dis=dist, k=5)
  testthat::expect_no_error(rg.test(E = E, n1 = num1, n2 = num2, k = 5, weigh.fun = weiMax, perm.num = 50))
  
  })

