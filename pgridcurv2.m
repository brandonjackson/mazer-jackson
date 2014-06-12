function [  ] = pgridcurv2( pf, offset, winsize, varargin )
%function [ ] = pgridcurv2(pf, offset, winsize, varargin)
%
% PGRIDCURV2 process and plots gridcurv data (hexagonal packed gabors)
%
%    [] = pgridcurv2(pf, offset, winsize, [options])
%
%  INPUT
%    file      - p2m PF struct
%    offset    - bin start time, in ms
%    winsize   - bin width, in ms
%
%    options - specified as <'option', value> pairs:
%      title     - title of plot
%      newfig    - (0/1) new figure window?
%      verbose   - (0/1) verbose plot?
%
%  OUTPUT
%    (figure)
%
%
% Wed Sep 18 10:00:00 2013 brandon
%   Prepped for deployment
%
% Mon Sep 23 11:05:26 2013 mazer 
%   - note I renamed this to pgridcurv2.m to avoid conflict with
%     the existing pgridcurv (which I think I wrote..)
%

p = inputParser;
addRequired(p,'pf');
addRequired(p,'offset');
addRequired(p,'winsize');
addParamValue(p,'title',[]);
addParamValue(p,'newfig',1);
addParamValue(p,'verbose',0);
try
  parse(p,pf,offset,winsize,varargin{:});
catch
  error('usage: pgridcurv2(pf, offset, winsize, ...opts..');
end

% optional params
if isempty(p.Results.title)
    figure_title = 'pgridcurv2';
else
    figure_title = p.Results.title;
end


GC = GridCurv(pf,offset,winsize);

orilist = [0,22,45,68,90,112,135,158]; % @TODO set dynamically

% Calculate Theta Mean
p_theta_mean = GC.orientation_means;
durations_dist = GC.durations_dist;

% Spike-triggered stimuli ensemble
sts_100 = GC.sts;

bootstats_100 = GC.generateStats();
bootstats_50 = GC.generateStats('fraction',0.5);
bootstats_nori = GC.generateStats('fraction',1/length(orilist));

% Calculate how standard error has changed over time
se_x = [(1/length(orilist))*100 50 100];
se_y = [std(bootstats_nori.raw(:)),...
        std(bootstats_50.raw(:)),...
        std(bootstats_100.raw(:))];
    
%% Conditionals

% Find cell and orientation with the max firing rate
max_count = max(bootstats_100.p_theta(:));
[max_cell,max_ori] = find(bootstats_100.p_theta==max_count);

[stats_max, stats_max_plus] = GC.generateConditionalStats(max_cell,max_ori);

[stats_center_0, stats_center_0_plus] = GC.generateConditionalStats('center',1);

if p.Results.verbose == 1
    [stats_center_45, stats_center_45_plus] = GC.generateConditionalStats('center',3);
    [stats_center_90, stats_center_90_plus] = GC.generateConditionalStats('center',5);
    [stats_center_135, stats_center_135_plus] = GC.generateConditionalStats('center',7);
end

%% Plots
if p.Results.newfig == 1
    f = figure('Name', figure_title, ...
              'NumberTitle', 'off', ...
              'Toolbar', 'none');
else
    f = gcf;
end


%% Verbose Plot
if p.Results.verbose == 1
    p = uiextras.TabPanel('Parent',f,'Padding',5);


    % Bars
    bars = uiextras.TabPanel('Parent',p,'Padding',5);

    bars_full = uipanel('Parent',bars,'Title','Full STS');
    GC.plotBars('Parent',bars_full,'stats',bootstats_100);

    bars_max = uipanel('Parent',bars,'Title','STS | maxori');
    subplot(1,2,1,'Parent',bars_max);
    GC.plotBars('stats',stats_max);
    subplot(1,2,2,'Parent',bars_max);
    GC.plotBars('stats',stats_max_plus);

    bars_center_0 = uipanel('Parent',bars,'Title','Center @ 0');
    subplot(1,2,1,'Parent',bars_center_0);
    GC.plotBars('stats',stats_center_0);
    subplot(1,2,2,'Parent',bars_center_0);
    GC.plotBars('stats',stats_center_0_plus);

    bars_center_1 = uipanel('Parent',bars,'Title','Center @ 45');
    subplot(1,2,1,'Parent',bars_center_1);
    GC.plotBars('stats',stats_center_45);
    subplot(1,2,2,'Parent',bars_center_1);
    GC.plotBars('stats',stats_center_45_plus);


    bars_center_2 = uipanel('Parent',bars,'Title','Center @ 90');
    subplot(1,2,1,'Parent',bars_center_2);
    GC.plotBars('stats',stats_center_90);
    subplot(1,2,2,'Parent',bars_center_2);
    GC.plotBars('stats',stats_center_90_plus);

    bars_center_3 = uipanel('Parent',bars,'Title','Center @ 135');
    subplot(1,2,1,'Parent',bars_center_3);
    GC.plotBars('stats',stats_center_135);
    subplot(1,2,2,'Parent',bars_center_3);
    GC.plotBars('stats',stats_center_135_plus);


    bars.TabNames = {'Full','Max','C @ 0', 'C @ 45', 'C @ 90', 'C @135'};
    bars.SelectedChild = 1;

    % STC Analysis

    real_stc = GC.stcAnalysis(GC.sts,GC.stimulus_avg);

    stc = uipanel('Parent',p,'title','STC Analysis');
    subplot(2,3,1,'Parent',stc);
    GC.plotVectorGabors(real_stc.sta);
    title('STA');

    subplot(2,3,2,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(1,:));
    title('Eigen Vector #1');

    subplot(2,3,3,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(2,:));
    title('Eigen Vector #2');

    if size(real_stc.eigen_vectors,1)>=3
        subplot(2,3,4,'Parent',stc);
        GC.plotVectorGabors(real_stc.eigen_vectors(3,:));
        title('Eigen Vector #3');
    end

    if size(real_stc.eigen_vectors,1)>=4
        subplot(2,3,5,'Parent',stc);
        GC.plotVectorGabors(real_stc.eigen_vectors(4,:));
        title('Eigen Vector #4');
    end

    subplot(2,3,6,'Parent',stc);
    scaled_real_eigenvalues = 100 * (real_stc.eigen_values / sum(real_stc.eigen_values(:)));

    [control_eigenvalues_mean, control_eigenvalues_std] = GC.findControlEigenvalues();
    scaled_control_eigenvalues = 100 * (control_eigenvalues_mean / sum(control_eigenvalues_mean(:)));
    scaled_control_eigenvalues_std = 100 * (control_eigenvalues_std / sum(control_eigenvalues_mean(:)));

    plot(scaled_real_eigenvalues,'-bo');
    hold on;
    errorbar(1:length(scaled_real_eigenvalues),scaled_control_eigenvalues,scaled_control_eigenvalues_std,'-r');
    legend('Real','Control','Location','NorthEast');

    %(scaled_control_eigenvalues,'--ro');
    ylabel('% variance explained');
    title('eigenvalues');


    % Tuning Curves
    curves = uiextras.TabPanel('Parent',p,'Padding',5);

    curves_full = uipanel('Parent',curves,'Title','Full Dataset');
    GC.plotTuningCurves('Parent',curves_full,'stats',bootstats_100);

    curves_half = uipanel('Parent',curves,'Title','Half of Dataset');
    GC.plotTuningCurves('Parent',curves_half,'stats',bootstats_50);

    curves_nori = uipanel('Parent',curves,'Title','1 / n_ori');
    GC.plotTuningCurves('Parent',curves_nori,'stats',bootstats_nori);

    curves.TabNames = {'Full','Half','1/(n_ori)'};
    curves.SelectedChild = 1;

    % Polar Plots
%     polars = uiextras.TabPanel('Parent',p,'Padding',5);

    polars_full = uipanel('Parent',p,'Title','Full Dataset');
    GC.plotPolars('Parent',polars_full,'stats',bootstats_100);
%     % plot_gridcurv_sts('polar',bootstats_100,'Parent',polars_full,'grid',GC.grid_xy);
% 
%     polars_half = uipanel('Parent',polars,'Title','Half of Dataset');
%     GC.plotPolars('Parent',polars_half,'stats',bootstats_50);
% 
%     polars_nori = uipanel('Parent',polars,'Title','1 / n_ori');
%     GC.plotPolars('Parent',polars_nori,'stats',bootstats_nori);
% 
% 
%     polars.TabNames = {'Full','Half','1/(n_ori)'};
%     polars.SelectedChild = 1;


    % Stimuli Statistics
    stimstats = uipanel('Parent',p,'Title','Stimulus Statistics');

    subplot(2,2,1,'Parent',stimstats);
    GC.plotCrossValidation();
    % plot(se_x, se_y,'-ro');
    % title('Standard Error');
    % xlabel('% of Dataset');
    % ylabel('Standard Error');
    % axis([0 100 0 max(se_y)*1.1]);

    subplot(2,2,2,'Parent',stimstats);
    
    % collapse across all cells to measure orientation anisotropy
    if size(p_theta_mean,1) > 1
        theta_mean_collapsed = mean(p_theta_mean,1);
        theta_std_collapsed = std(p_theta_mean,1);
    % but don't collapse if there was only one cell!
    else
        theta_mean_collapsed = p_theta_mean;
        theta_std_collapsed = zeros(size(p_theta_mean));
    end

    errorbar(orilist,theta_mean_collapsed,theta_std_collapsed);
    set(gca,'XTick',orilist);
    hold on;
    plot(-20:180,ones(201,1)*0.125,'-r');
    ylabel('p(\theta)');
    axis([-20 180 0.0625 0.1875]);
    title('Time per Orientation');

    subplot(2,2,3,'Parent',stimstats);
    hist(durations_dist);
    xlabel('ms');
    ylabel('n');
    title(['Durations (total=' num2str(sum(durations_dist)) 'ms, mean=' num2str(mean(durations_dist)) 'ms)']);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 .63 1],'EdgeColor','w');

    subplot(2,2,4,'Parent',stimstats);
    GC.plotRF();

%     % Spike Statistics
%     spikestats = uipanel('Parent',p,'Title','Spike Statistics');
%     psth = axes('Parent',spikestats);
%     phist(GC.pf);

    p.TabSize = 80;
    p.TabNames = { 'Bar Plots','STC Analysis', 'Tuning Curves','Polar Plots', 'Stimuli Stats'  };
    p.SelectedChild = 1;
    
%     h = figure();
%     
%     p1 = uipanel('parent',h,'Position',[0 0 0.5 1],'Title','Full STA');
%     GC.plotPolars('Parent',p1,'stats',bootstats_100);
% 
%     p2 = uipanel('parent',h,'Position',[0.5 0 0.5 1],'Title','STA | Optimal Orientation');
%     GC.plotPolars('parent',p2,'stats',stats_center_0);
%     
% 
%     h2 = figure();
%     
%     p3 = uipanel('parent',h2,'Position',[0 0 0.5 1],'Title','Full STA');
%     GC.plotBars('Parent',p3,'stats',bootstats_100);
% 
%     p4 = uipanel('parent',h2,'Position',[0.5 0 0.5 1],'Title','STA | Optimal Orientation');
%     GC.plotBars('parent',p4,'stats',stats_center_0);

    
%% Simpler Plot
% this plot uses a single level of tab panels, and is here for posterity
% in case anyone ever wants to use it.
elseif p.Results.verbose == 2
    
    p = uiextras.TabPanel('Parent',f,'Padding',5);


    % Bars
    bars = uipanel('Parent',p,'Title','Full STS & Center @ 0 +/-');
    
    subplot(1,2,1,'Parent',bars);
    GC.plotBars('stats',bootstats_100);
    colorbar();
    subplot(1,2,2,'Parent',bars);
    GC.plotBars('stats',stats_center_0_plus);
    colorbar();

    % STC Analysis

    real_stc = GC.stcAnalysis(GC.sts,GC.stimulus_avg);

    stc = uipanel('Parent',p,'title','STC Analysis');
    subplot(2,3,1,'Parent',stc);
    GC.plotVectorGabors(real_stc.sta);
    title('STA');

    subplot(2,3,2,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(1,:));
    title('Eigen Vector #1');

    subplot(2,3,3,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(2,:));
    title('Eigen Vector #2');

    subplot(2,3,4,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(3,:));
    title('Eigen Vector #3');

    subplot(2,3,5,'Parent',stc);
    GC.plotVectorGabors(real_stc.eigen_vectors(4,:));
    title('Eigen Vector #4');

    subplot(2,3,6,'Parent',stc);
    scaled_real_eigenvalues = 100 * (real_stc.eigen_values / sum(real_stc.eigen_values(:)));

    [control_eigenvalues_mean, control_eigenvalues_std] = GC.findControlEigenvalues();
    scaled_control_eigenvalues = 100 * (control_eigenvalues_mean / sum(control_eigenvalues_mean(:)));
    scaled_control_eigenvalues_std = 100 * (control_eigenvalues_std / sum(control_eigenvalues_mean(:)));

    plot(scaled_real_eigenvalues,'-bo');
    hold on;
    errorbar(1:length(scaled_real_eigenvalues),scaled_control_eigenvalues,scaled_control_eigenvalues_std,'-r');
    legend('Real','Control','Location','NorthEast');

    %(scaled_control_eigenvalues,'--ro');
    ylabel('% variance explained');
    title('eigenvalues');


    % Tuning Curves

    curves = uipanel('Parent',p,'Title','Tuning Curves (Full Dataset)');
    GC.plotTuningCurves('Parent',curves,'stats',bootstats_100);

    % Polar Plots

    polars = uipanel('Parent',p,'Title','Polar Plots (Full Dataset)');
    GC.plotPolars('Parent',polars,'stats',bootstats_100);

    % Stimuli Statistics
    stimstats = uipanel('Parent',p,'Title','Stimulus Statistics');

    subplot(2,2,1,'Parent',stimstats);
    GC.plotCrossValidation();

    subplot(2,2,2,'Parent',stimstats);
    theta_mean_collapsed = mean(p_theta_mean,1);
    theta_std_collapsed = std(p_theta_mean,1);
    errorbar(orilist,theta_mean_collapsed,theta_std_collapsed);
    set(gca,'XTick',orilist);
    hold on;
    plot(-20:180,ones(201,1)*0.125,'-r');
    ylabel('p(\theta)');
    axis([-20 180 0.0625 0.1875]);
    title('Time per Orientation');

    subplot(2,2,3,'Parent',stimstats);
    hist(durations_dist);
    xlabel('ms');
    ylabel('n');
    title(['Durations (total=' num2str(sum(durations_dist)) 'ms, mean=' num2str(mean(durations_dist)) 'ms)']);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 .63 1],'EdgeColor','w');

    subplot(2,2,4,'Parent',stimstats);
    GC.plotRF();



    p.TabSize = 100;
    p.TabNames = { 'Bar Plots','STC Analysis', 'Tuning Curves','Polar Plots', 'Stimuli Stats'};
    p.SelectedChild = 1;

%% Simplest Plot
else
     % Tuning Curves
    curves = uipanel('Parent',f,'Title',['pgridcurv2  ' num2str(offset) '-' num2str(offset+winsize) 'ms  ' GC.pf.src],'Position',[0.025 0.42 .95 0.55]);
    GC.plotTuningCurves('Parent',curves,'stats',bootstats_100);
    
     % Bars
    bars = uipanel('Parent',f,'Title','Full STS, Center @ 0 +/-, and Cross-Validation','Position',[0.025 0.025 0.95 0.37]);
    
    a1 = axes('Parent',bars,'Position',[0 0 0.34 1]);
    GC.plotBars('stats',bootstats_100,'Parent',a1);
    colorbar('location','WestOutside');
    a2 = axes('Parent',bars,'Position',[0.34 0 0.34 1]);
    GC.plotBars('stats',stats_center_0_plus,'Parent',a2);
    colorbar('location','WestOutside');
    a3 = axes('Parent',bars,'Position',[0.72 0.15 0.26 0.75]);
    GC.plotCrossValidation('Parent',a3);
    ylabel('');
end


