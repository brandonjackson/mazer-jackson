% pgridcurv2
% This function is the interface to the GridCurv class and displays a
% dashboard-like summary for each cell, a la pgratrev.
% It is called pgridcurv2 since Jamie wrote pgridcurv before I pushed my
% changes to the svn repo :-)

% Panel A: Simple Plot
% The default plot is simple, and doesn't use tabs to make it
% linux-friendly
pgridcurv2(pffind('romeo0284*gridcurv*002'),50,50);

% Panel B: Verbose
pgridcurv2(pffind('romeo0284*gridcurv*002'),50,50,'verbose',1);
