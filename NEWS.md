Changes in version 0.1.1

Improvements

    - We added a new parameter called 'return_model' to the cv.FastStepGraph function to specify when to return the fitted model, i.e. internally calls FastStepGraph with the optimal alpha_f and alpha_b parameters and returns Omega.

    - We added a new parameter called b_coef with initial value 0.5. This parameter serves to control the empirical rule alpha_b=b_coef*alpha_f during the initial search for the optimal alpha_f parameter while alpha_b remains fixed, after finding optimal alpha_f, alpha_b is varied to find its ideal value. The default value of b_coef is 0.5, but with real data sets it may be necessary to use a lower value.

Bug fixes

    - In the previous version only one value of alpha_b was tested, now several different values of this parameter are tested to find the ideal value that produces the minimum loss during cross-validation.

# FastStepGraph 0.1.0

* Initial CRAN submission.
