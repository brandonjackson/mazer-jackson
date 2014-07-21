% Find Kernel Correlation / Performance
% How strong is the correlation between the firing rates predicted by the
% kernel and the observed rates?

AP = AnglePlay(pffind('romeo0300*curvplay'),60,100);

% Panel A: Full Bootstrapped Analysis
% Uses 90% of data as training set, resampled 100 times.
AP.findKernelCorrelation(0.9,'bootstrap',1);

% Panel B: Single Analysis (no bootstrap)
% Uses 100% of data as training set, generates scatter plot.
AP.findKernelCorrelation(1);