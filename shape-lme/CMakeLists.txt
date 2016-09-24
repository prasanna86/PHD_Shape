## shape-lme binaries
ADD_EXECUTABLE(hypothesis-testing-tutorial src/hypothesis-testing-tutorial.cxx)
TARGET_LINK_LIBRARIES(hypothesis-testing-tutorial ${ARMADILLO_LIBRARIES})

ADD_EXECUTABLE(estimate-covariate-mixed-effects src/estimate-covariate-mixed-effects.cxx) 
TARGET_LINK_LIBRARIES(estimate-covariate-mixed-effects ${ARMADILLO_LIBRARIES})
