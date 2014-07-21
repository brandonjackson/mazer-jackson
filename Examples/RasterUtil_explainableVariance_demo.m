% Explainable Variance

% PANEL A
% If jackknife arg not provided, then it does a jackknife by default.
% It then displays histogram of jackknifed explainable r's.
RasterUtil.explainableVariance(pffind('bert0270*curvplay'),80,50);

% PANEL B
% If jackknife set to false, then it displays a scatter plot
RasterUtil.explainableVariance(pffind('bert0270*curvplay'),80,50,'jackknife',0);