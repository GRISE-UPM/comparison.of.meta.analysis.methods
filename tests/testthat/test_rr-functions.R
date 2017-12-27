# ----------------------------------------------------------------------------------------------------------------
# Research             : A procedure for selecting an appropriate meta-analysis method: A comparison of approaches
# File name            : test_rr-functions.R
# File type            : Unit tests
# Purpose              : Perform the unit testing of the functions included in the file rr-functions.R
# Creation date        : History available at GitHub
# R Script contributors: Oscar Dieste, Omar S. Gomez
# ----------------------------------------------------------------------------------------------------------------
#

context('Response ratio tests')

# Unit tests for: tell_rr
#

test_that('the rr of control = 1, treatment= 1, not using logs is 1', {
     CONTROL <- 1
     TREATMENT <- 1
     USELOGS <- FALSE
     RESULT <- 1 # TREATMENT / CONTROL
     expect_equal(RESULT, tell_rr(CONTROL, TREATMENT, USELOGS))
})

test_that('the rr of control = 1, treatment= 1, when using logs is 0', {
     CONTROL <- 1
     TREATMENT <- 1
     USELOGS <- TRUE
     RESULT <- 0 # log(TREATMENT / CONTROL)
     expect_equal(RESULT, tell_rr(CONTROL, TREATMENT, USELOGS))
})

test_that('the variance of a rr is ...incomplete test', {
     fail()
})
