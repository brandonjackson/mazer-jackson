% Compare Gratrev and Angleplay Model Performance

% AnglePlay Stimulus
ap_pf = pffind('romeo0295*curvplay');
ap_offset = 50;
ap_winsize = 100;

% GratRev Grating Stimulus
gr_pf = pffind('romeo0295*gratrev*001');
gr_offset = 70;
gr_winsize = 100;

% Compare Models
[gr_gr_score,gr_ap_score,ap_ap_score] = AnglePlay.gratrevModelComparison(ap_pf,...
    ap_offset,ap_winsize,gr_pf,gr_offset,gr_winsize);
