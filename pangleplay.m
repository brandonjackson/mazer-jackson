function [  ] = pangleplay( pf, latency, winsize, varargin )
%function [ ] = pangleplay(pf, latency, winsize, varargin)
%
% PANGLEPLAY process and plots parabola/abs curvplay data
%
%    [] = pangleplay(pf, latency, winsize, [options])
%
%  INPUT
%    file      - p2m PF struct
%    latency   - temporal latency (ms)
%    winsize   - temporal integration window (ms)
%
%    options - specified as <'option', value> pairs:
%      significant - (0/1) show only significant values (default=1)
%      title     - title of plot
%      newfig    - (0/1) new figure window?
%      imagedir  - (string) path to stimulus images
%
%  OUTPUT
%    (figure)
%
%
% Wed Sep 18 11:30:00 2013 brandon
%   building interface for deployment

p = inputParser;
addRequired(p,'pf');
addRequired(p,'latency');
addParamValue(p,'title',[]);
addParamValue(p,'imagedir',0);
addParamValue(p,'significant',1);
addParamValue(p,'newfig',1);
addParamValue(p,'printable',0);
addParamValue(p,'max_rate',0);
parse(p,pf,latency,varargin{:});

AP = AnglePlay(pf,latency,winsize);

% optional params
if isempty(p.Results.title)
    figure_title = ['pangleplay ' AP.exper ' ' num2str(latency) '-' num2str(latency+winsize) 'ms'];
else
    figure_title = p.Results.title;
end


if p.Results.newfig == 1
    f = figure('Name', figure_title, ...
              'NumberTitle', 'off', ...
              'Toolbar', 'none');
else
    f = gcf;
end

subplot(3,3,1);
SuperUtil.plotPsthWindow(AP.pf, latency, winsize);
phist(AP.pf,'pre',0,'post',300,'suptitle',0);

subplot(3,3,2);
exper_parts = strsplit(AP.exper,'.');
cell_name = exper_parts{1};
% gratrev_str = [exper_parts{1} '*gratrev'];

try
    results_filename = ['/lab/data/gratrev_results/' cell_name '_gratrev_results'];
    load(results_filename);
    tuning = gratrev_results.tuning;
    errorbar(tuning.x,tuning.y,tuning.e,'-k');
    xlim([0 max(tuning.x)]);
    ylim([0 (max(tuning.y) + 0.333*max(tuning.y))]);
    ylabel('s/s');
    xlabel('ori');
    title(['Gratrev Tuning (lat=' num2str(gratrev_results.latency) ', winsize=' num2str(gratrev_results.winsize) ')'],'fontweight','bold');
catch err
    title('Gratrev Tuning Unavailable','fontweight','bold');
end

subplot(3,3,3);
AP.plotCrossValidation();

subplot(3,3,4);
AP.plotAbsTuning('unit','angles','zscore',0,'significant',p.Results.significant);
title('Abs Tuning','fontweight','bold');
abs_rates = caxis;

subplot(3,3,5);
AP.plotParabolaTuning('unit','angles','zscore',0,'significant',p.Results.significant);
title('Parabola Tuning','fontweight','bold');
par_rates = caxis;

subplot(3,3,6);
AP.plotTypes();

subplot(3,3,7);
if p.Results.imagedir ~= 0
    AP.plotImageDomain('show_final',1,'plot_parabola',0,'plot_abs',1,'imagedir',p.Results.imagedir,'significant',p.Results.significant);
else
    AP.plotImageDomain('show_final',1,'plot_parabola',0,'plot_abs',1,'significant',p.Results.significant);
end
title('Abs STA (Contrast)','fontweight','bold');

subplot(3,3,8);
if p.Results.imagedir ~= 0
    AP.plotImageDomain('show_final',1,'plot_parabola',1,'plot_abs',0,'imagedir',p.Results.imagedir,'significant',p.Results.significant);
else
    AP.plotImageDomain('show_final',1,'plot_parabola',1,'plot_abs',0,'significant',p.Results.significant);
end
title('Parabola STA (Contrast)','fontweight','bold');

subplot(3,3,9);
if p.Results.imagedir ~= 0
    AP.plotImageDomain('show_final',1,'imagedir',p.Results.imagedir,'significant',p.Results.significant);
else
    AP.plotImageDomain('show_final',1,'significant',p.Results.significant);
end
title('Full STA (Contrast)','fontweight','bold');


% if max_rate set, scale all tuning plots to have the same color scale
% this is used when doing batch pangleplay plots, comparing across multiple
% time scales
if p.Results.max_rate ~= 0
    subplot(3,3,4);
    caxis([0 p.Results.max_rate]);
    subplot(3,3,5);
    caxis([0 p.Results.max_rate]);
    subplot(3,3,6);
    line([0 p.Results.max_rate],[0 p.Results.max_rate],'Color','k');
    xlim([0 p.Results.max_rate]);
    ylim([0 p.Results.max_rate]);
else
    % set firing rate color scales to be the same for both parabola and
    % absolute value tuning plots

    if(abs_rates(2) > par_rates(2))
        subplot(3,3,4);
        caxis(abs_rates);
        subplot(3,3,5);
        caxis(abs_rates);
    else
        subplot(3,3,4);
        caxis(par_rates);
        subplot(3,3,5);
        caxis(par_rates);
    end
end

boxtitle(['pangleplay   ' num2str(latency) '-' num2str(latency+winsize) 'ms   ' AP.pf.src]);

end


