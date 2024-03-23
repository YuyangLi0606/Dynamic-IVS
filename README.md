# Dynamic-IVS
The coding supplemental material for paper "Modelling Implied Volatility Surface Using B-splines with Time-dependent Coefficients Predicted by Tree-based Machine Learning Methods".

The algorithm could be divided into three steps:

- Constructing bivariate B-spline basis (in R): construct B-spline basis with the help with the help of *mgcv* package. Fit the IVS on a daily scheme and estimate $\alpha_t$'s.
- Fitting the time series of the dynamic coefficients: fit a Gradient Boosted Regression Tree and predict the future $\alpha$'s with the historical information.
- Forecasting future IVS: plug the predicted $\alpha$ into the common (fixed) basis and obtain the forecasted IVS.
