## shape-lme binaries
ADD_EXECUTABLE(hypothesis-testing-tutorial src/hypothesis-testing-tutorial.cxx)
TARGET_LINK_LIBRARIES(hypothesis-testing-tutorial ${ARMADILLO_LIBRARIES})

ADD_EXECUTABLE(mixed-effects-covariate-estimation src/mixed-effects-covariate-estimation.cxx) 
TARGET_LINK_LIBRARIES(mixed-effects-covariate-estimation ${ARMADILLO_LIBRARIES})
