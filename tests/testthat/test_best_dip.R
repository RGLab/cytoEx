test_that("Test best.dip singlet selection unchanged: ", {
    set.seed(172384)


    res <- best.dip(fr.singlet,parallel_type="multicore",mc.cores=4)

    #this data set should have two channels with non-zero p-value for dip-statistic
    expect_true(all(res[c(2,6),intial.p.values] > 4.5e-3))
    #data set should select "<B710-A>"
    expect_true(res[channel=="<B710-A>", area.ratio] == 1)
})



