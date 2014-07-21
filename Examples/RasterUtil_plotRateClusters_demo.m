% Plot Rate Clusters

% Uses k-means clustering to look for clusters in the firing rates to
% different stimuli. By default it scales/normalizes the rates so that the
% highest firing rate for each is equal to one.
n_clusters = 2;
transient_onset_time = 50;
RasterUtil.plotRateClusters(pffind('romeo0300*gratrev*003'),...
    n_clusters,...
    transient_onset_time);