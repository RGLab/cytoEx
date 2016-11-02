test_that("Test best.dip.ICL across several data sets. ", {
    set.seed(1239214)
    
    #
    #tests for singlets
    #

    res <- best.dip.ICL(fr.singlet,parallel_type="multicore",mc.cores=4)

    #this data set should have two channels with non-zero p-value for dip-statistic
    expect_true(all(res[c(2,6),intial.p.values] > 4.5e-3))
    expect_equal(sum(res[,"intial.p.values",with=FALSE] > 0),2)
   
    #the initial p-value of the dip statistic for two channels should not change with seed set.
    expect_equal(round(res[channel=="<R660-A>",intial.p.values],3), 0.005)
    expect_equal(round(res[channel=="<G560-A>",intial.p.values],3), 1.000)
    
    #this method should select "<V450-A>"
    expect_true(res[channel=="<V450-A>", area.ratio] == 1)

    #the method should then reject all other channels (set to zero)
    expect_equal(res[!(channel=="<V450-A>"), area.ratio],c(0,0,0,0,0,0))
    
    #it passes over <B710-A> because of negative difference in ICL between 2 and 3 mode
    expect_true(res[channel=="<B710-A>", diff.icl] < 0)

    #it doesn't compute difference in ICL for <G560-A>
    expect_true(res[channel=="<G560-A>", diff.icl]== -Inf)

    #difference in icl should be computed for the rest
    expect_equal(min(sapply(res[!(channel == "<G560-A>"),diff.icl],is.finite)),1)

    #it ranks seven channels
    expect_equal(res[, score], c(4,2,6,7,3,1,5))

    #the sub-sampled dip-statistic p-value for  <V545-A> should exceed 0.05.
    expect_true(res[channel=="<V545-A>",second.p.values] > 0.05)

    #the significance level (with bonferonni correction) shouldn't change
    expect_equal(max(round(res[,"b.alpha",with=FALSE],3)),0.007)
    expect_equal(min(round(res[,"b.alpha",with=FALSE],3)),0.007)


    #
    #tests for singlets/CD3+
    #

    res.cd3 <- best.dip.ICL(fr.cd3,parallel_type="multicore",mc.cores=4)

    #results should have four channels with non-zero p-value for dip-statistic
    expect_true(all(res.cd3[c(2,4,5,6),intial.p.values] > 0.01))
    expect_equal(sum(res.cd3[,"intial.p.values",with=FALSE] > 0),4)
   
    #the initial p-value of the dip statistic for two channels should not change with seed set.
    expect_equal(round(res.cd3[channel %in% c("<R660-A>","<V450-A>","<V545-A>","<G560-A>"),intial.p.values],3),
                 c(1.000,1.000,1.000,0.738))
    
    #method should select "<B710-A>"
    expect_true(res.cd3[channel=="<B710-A>", area.ratio] == 1)

    #the method should then reject all other channels (set to zero)
    expect_equal(res.cd3[!(channel=="<B710-A>"), area.ratio],c(0,0,0,0,0,0))
    

    #it shouldn't compute difference in ICL for the following
    expect_true(max(res.cd3[channel %in% c("<R660-A>","<V450-A>","<V545-A>","<G560-A>"), diff.icl]) == -Inf)

    #difference in icl should be computed for the rest
    expect_equal(min(sapply(res.cd3[!(channel %in% c("<R660-A>","<V450-A>","<V545-A>","<G560-A>")),diff.icl],
                            is.finite)),1)

    
    #it ranks seven channels
    expect_equal(res.cd3[, score], c(7,1,6,2,3,4,5))

    #the sub-sampled dip-statistic p-value should all be below 0.01
    expect_true(max(res.cd3[channel %in% c("<B710-A>","<R780-A>","<G780-A>"),second.p.values]) < 0.01)

    #the significance level (with bonferonni correction) shouldn't change
    expect_equal(max(round(res.cd3[,"b.alpha",with=FALSE],3)),0.007)
    expect_equal(min(round(res.cd3[,"b.alpha",with=FALSE],3)),0.007)

    #
    #tests for singlets/CD3+/CD8
    #

    res.cd8 <- best.dip.ICL(fr.cd8,parallel_type="multicore",mc.cores=4)

    #results should have five channels with non-zero p-value for dip-statistic
    expect_true(all(res.cd8[c(1,2,3,4,5),intial.p.values] > 0.1))
    expect_equal(sum(res.cd8[,"intial.p.values",with=FALSE] > 0),5)
   
    #the initial p-value of the dip statistic for five channels should not change with seed set.
    expect_equal(round(res.cd8[channel %in% c("<B710-A>","<R660-A>","<R780-A>","<V450-A>","<V545-A>"),intial.p.values],3),
                 c(0.995,0.998,1.000,1.000,0.996))


    #Shouldn't compute ICL for any of these
    expect_true(max(res.cd8[channel %in% c("<B710-A>","<R660-A>","<R780-A>","<V450-A>","<V545-A>"),diff.icl]) == -Inf)
    
    #difference in icl should be computed for the rest
    expect_equal(min(sapply(res.cd8[!(channel %in% c("<B710-A>","<R660-A>","<R780-A>","<V450-A>","<V545-A>")),diff.icl],
                            is.finite)),1)

    
    
    #method should select "<G560-A>"
    expect_true(res.cd8[channel=="<G560-A>", area.ratio] == 1)

    #the method should then reject all other channels (set to zero)
    expect_equal(res.cd8[!(channel=="<G560-A>"), area.ratio],c(0,0,0,0,0,0))
    

    #it ranks seven channels
    expect_equal(res.cd8[, score], c(1,2,3,4,5,7,6))

    #the sub-sampled dip-statistic p-value should be small for the winning channel, large for the losing
    #since precise value will change if the seed changes, perform appoximate check
    expect_true(res.cd8[channel == "<G560-A>",second.p.values] < 0.2)
    expect_true(res.cd8[channel == "<G780-A>",second.p.values] > 0.2)

    #once again, the significance level (with bonferonni correction) shouldn't change
    expect_equal(max(round(res.cd8[,"b.alpha",with=FALSE],3)),0.007)
    expect_equal(min(round(res.cd8[,"b.alpha",with=FALSE],3)),0.007)
    
    
    #
    #tests for DC::Monocytes/Live
    #

    res.dc.mon <- best.dip.ICL(fr.dc.mon,parallel_type="multicore",mc.cores=4)

    #it ranks eight channels
    expect_equal(res.dc.mon[, score], c(3,4,1,7,2,5,8,6))

    #method should select "CD14"
    expect_true(res.dc.mon[channel=="CD14", area.ratio] == 1)

    #the method should then reject all other channels (set to zero)
    expect_equal(res.dc.mon[!(channel=="CD14"), area.ratio],c(0,0,0,0,0,0,0))

    #the initial p-value of the dip statistic for two channels should not change across tests
    expect_equal(round(res.dc.mon[channel %in% c("CD123","CD16"),intial.p.values],3),c(0.997,0.008))

    #the difference in icl should not be computed for the two channels
    expect_equal(res.dc.mon[channel %in% c("CD123","CD16"),diff.icl],c(-Inf,-Inf))

    #difference in icl should be computed for the rest
    expect_equal(min(sapply(res.dc.mon[!(channel %in% c("CD123","CD16")),diff.icl],is.finite)),1)
    
    #CD56 and Lineage should have large second.p.values
    #since the exact value will change if set.seed is not invoked, perform an approximate check.
    expect_true(res.dc.mon[channel=="CD56",second.p.values] > 0.25)
    expect_true(res.dc.mon[channel=="Lineage",second.p.values] > 0.25)


    #significance level should be the same for all rows
    expect_equal(max(round(res.dc.mon[,"b.alpha",with=FALSE],3)),0.006)
    expect_equal(min(round(res.dc.mon[,"b.alpha",with=FALSE],3)),0.006)


    #
    #tests for bcell Live
    #
    res.bcell <- best.dip.ICL(fr.bcell,parallel_type="multicore",mc.cores=4)

    #it ranks eight channels
    expect_equal(res.bcell[, score], c(1,7,8,6,2,5,3,4))

    #method should select "CD19"
    expect_true(res.bcell[channel=="CD19", area.ratio] == 1)

    #the method should then reject all other channels (set to zero)
    expect_equal(res.bcell[!(channel=="CD19"), area.ratio],c(0,0,0,0,0,0,0))
        
    #the initial p-value of the dip statistic for Live channel should not change across tests
    expect_equal(round(res.bcell[channel == "Live",intial.p.values],3),c(0.998))

    #difference in icl shouldn't be computed for Live
    expect_equal(res.bcell[channel == "Live",diff.icl],c(-Inf))

    #difference in icl should be computed for the rest
    expect_equal(min(sapply(res.bcell[!(channel=="Live"),diff.icl],is.finite)),1)

    #difference in icl should be negative for IgD
    expect_true(res.bcell[channel == "IgD",diff.icl]<0)
    

    #significance level should be the same for all rows
    expect_equal(max(round(res.bcell[,"b.alpha",with=FALSE],3)),0.006)
    expect_equal(min(round(res.bcell[,"b.alpha",with=FALSE],3)),0.006)


    #IgD and CD24 should have large second.p.values
    #since the exact value will change if set.seed is not invoked, perform an approximate check.
    expect_true(res.bcell[channel=="IgD",second.p.values] > 0.25)
    expect_true(res.bcell[channel=="CD24",second.p.values] > 0.25)

    }
)
