context("test exaustive gating..")

test_that("bcell ", {
  gs2 <- clone(gs.bcell, isNew = FALSE, isEmpty = FALSE)
  gs2 <- flowIncubator::swapChannelMarker(gs2)
  for(node in getChildren(gs2[[1]], "lymph"))
    Rm(node, gs2)

  min.percent <- 0.2
  openCyto::gating("lymph", gs2
                   , min.count = 2000
                   , min.percent = min.percent
                   , debug.mode = T
                   , parallel_type=parallel_type, mc.cores=mc.cores,  P.ITERS = P.ITERS
                   , gating.function = openCyto::mindensity
                   , marker.selection.function = best.dip.ICL
                   , marker.selection.args = list(P.ITERS=100)
  )

  expect_true(setequal(getNodes(gs2)[-(1:3)], c('/boundary/nonDebris/lymph/CD19+'
                                        , '/boundary/nonDebris/lymph/CD19-'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27+'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27+/CD38+'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27+/CD38-'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27+/CD38-/Live+'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27+/CD38-/Live-'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3+/CD27-'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3-'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3-/CD38+'
                                        , '/boundary/nonDebris/lymph/CD19-/CD3-/CD38-')
                  )
              )

})