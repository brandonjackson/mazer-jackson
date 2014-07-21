% Combined PSTH of Different Experiments for a Cell

% Only Plotting 5 tasks, to make plot more readable. Full list attached
% below for reference.
pf_list = cell(5,1);
pf_list{1} = pffind('romeo0300.spotmap.000');
pf_list{2} = pffind('romeo0300.gratrev.002');
pf_list{3} = pffind('romeo0300.gratrev.008');
pf_list{4} = pffind('romeo0300.gridcurv.006');
pf_list{5} = pffind('romeo0300.curvplay.005');

RasterUtil.plotExperPSTH(pf_list);

% % Full List
% romeo0300.curvplay.005
% romeo0300.gratrev.001
% romeo0300.gratrev.002
% romeo0300.gratrev.003
% romeo0300.gratrev.004
% romeo0300.gratrev.008
% romeo0300.gratrev.009
% romeo0300.gridcurv.006
% romeo0300.gridcurv.007
% romeo0300.spotmap.000