test_that("Test best.dip singlet selection unchanged: ", {
    set.seed(1239214)

    res <- best.dip.ICL(fr.singlet,parallel_type="multicore",mc.cores=4)

    #this data set should have two channels with non-zero p-value for dip-statistic
    expect_true(all(res[c(2,6),intial.p.values] > 4.5e-3))

    #this method should select "<V450-A>"
    expect_true(res[channel=="<V450-A>", area.ratio] == 1)

    #it passes over <B710-A> because of negative difference in ICL between 2 and 3 mode
    expect_true(res[channel=="<B710-A>", diff.icl] < 0)

    #it doesn't compute difference in ICL for <G560-A>
    expect_true(res[channel=="<G560-A>", diff.icl]== -Inf)

    #it ranks seven channels
    expect_equal(res[, score], c(4,2,6,7,3,1,5))
}
)
