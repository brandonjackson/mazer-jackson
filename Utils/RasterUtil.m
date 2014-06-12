classdef RasterUtil < SuperUtil
    % RASTERUTIL contains utilities related to spike trains and PSTHs
    % It is intended to complement / replace the functionality contained in
    % prast.
    
    properties
    end
    
    properties (Constant)
        % sigma of Guassian used in computing spike density functions
        SPIKE_DENSITY_SIGMA = 5
    end
    
    methods
        function RU = RasterUtil()
        end
    end
    
    methods(Static)
        
        function [sdf] = spikeDensity(spike_vector)
        % SPIKEDENSITY computes spike densities
        % This function converts the spike train SPIKE_VECTOR into a vector
        % of spike densities, SDF. SPIKE_VECTOR is assumed to be sampled at
        % 1ms resolution, and is often the raster from prast. The spike train 
        % is convolved with a gaussian with sigma taken from 
        % RasterUtil.SPIKE_DENSITY_SIGMA constant. The resulting vector 
        % is of the range [0,1], and is useful for calculating PSTHs
            sigma = RasterUtil.SPIKE_DENSITY_SIGMA; % in ms
            kernel_x = -3*sigma:3*sigma;
            kernel = normpdf(kernel_x,0,sigma);
            sdf = conv(spike_vector,kernel,'same');
        end
        
        function triggers = uniqueTriggers(pf)
            
            N_RECS = length(pf.rec);
            prefix = PFUtil.framePrefix(pf);
            
            triggers_cell = cell(N_RECS,1);
            for i=1:N_RECS
                [ix,~] = PFUtil.findEvents(pf,i,prefix);
                triggers_cell{i} = pf.rec(i).ev_e(ix);
                if size(triggers_cell{i},1) > size(triggers_cell{i},2)
                    triggers_cell{i} = triggers_cell{i}';
                end
            end
            triggers_cell_orig = triggers_cell;
            if size(triggers_cell{1},1) > size(triggers_cell{1},2)
                cat_dim = 1;
            else
                cat_dim = 2;
            end
            % concatenate all the nested cell arrays together
            while any(cellfun('isclass',triggers_cell,'cell'))
                triggers_cell = cat(cat_dim,triggers_cell{:});
            end
            
            triggers = unique(triggers_cell);
        end
        
        function [rates,ts,triggers] = rates(pf, pre, post)
        % RATES gets PSTHs (rates over time) for each unique stimulus frame
            
            if nargin < 3
                pre = -200;
                post = 400;
            end
            raster = prast(pf,'pre',pre,'post',post);
            ts = raster.time;
            
            [triggers, ~, trigger_is] = unique(raster.triggers);
            rates = zeros(length(triggers), size(raster.data,2));
            
            for i=1:length(triggers)
                mean_spike_vector = nanmean(raster.data(trigger_is==i,:),1);
                
                % Create spike density function
                sdf = RasterUtil.spikeDensity(mean_spike_vector);
                rates(i,:) = sdf * 1000;
            end
            max_rate = nanmean(rates(:)) + 4*nanstd(rates(:));
            
            if nargout < 1
                figure();
                for i=1:20
                    subplot(5,4,i);
                    area(ts, rates(randi(length(triggers)),:));
                    ylim([0 max_rate]);
                    xlim([0 250]);
                end
            end
        end
        
        function [explainable_r] = explainableVariance(pf,lat,winsize,varargin)
        % EXPLAINABLEVARIANCE measures correlation b/w odd and even trials
        % It takes the standard PF/LAT/WINSIZE arguments and returns the
        % pearson coefficient EXPLAINABLE_R. To find the explainable
        % variance this should be squared.
        %
        % If no output args specified then plots a scatter plot of the even
        % and odd trials, as well as a histogram of repetitions.
            
            p = inputParser;
            addRequired(p,'pf');
            addRequired(p,'lat');
            addRequired(p,'winsize');
            addParameter(p,'jackknife',1);
            parse(p,pf,lat,winsize,varargin{:});
            
            % Jackknife to estimate confidence intervals
            % (note: jackknifing used instead of bootstrap because the
            % repetitions used in bootstrapping introduce spurious
            % correlations)
            if p.Results.jackknife==1
                N_DRAWS = 10;
                [rs_mean,rs_err,rs_raw] = PFUtil.jackknife(pf,@RasterUtil.explainableVariance,N_DRAWS,lat,winsize,'jackknife',0);
                explainable_r = rs_mean;
                if nargout < 1
                    figure();
                    hist(rs_raw);
                    title(sprintf('Explainable Variance Distribution (mean r^2=%.2f)',rs_mean^2));
                    xlabel('r');
                end
                return;
            end
            
            raster = prast(pf);
            ts = raster.time;
            
            min_i = find(ts==lat);
            max_i = find(ts==(lat+winsize));
            
            [triggers, ~, trigger_is] = unique(raster.triggers);
            rates = zeros(length(triggers), size(raster.data,2));
            
            odds = [];
            evens = [];
            reps = zeros(length(triggers),1);
            
            for i=1:length(triggers)
                matches = raster.data(trigger_is==i,:);
                sdfs = zeros(size(matches));
                for j=1:size(matches,1)
                    sdfs(j,:) = RasterUtil.spikeDensity(matches(j,:));
                end
                
                % slice sdfs, leaving behind only the desired window of
                % data
                sdfs = sdfs(:,min_i:max_i);
                
                % convert to rates, then get mean rate
                rates = sdfs * 1000;
                mean_rates = mean(rates,2);
                
                % Computes mean of odd and even trials
                stim_odds = [];
                stim_evens = [];
                for j=1:10
                    if size(matches,1) >= (j*2)
%                         odds = [odds; matches(2*(j-1)+1)];
%                         evens = [evens; matches(2*(j-1)+2)];
                        stim_odds = [stim_odds; mean_rates(2*(j-1)+1)];
                        stim_evens = [stim_evens; mean_rates(2*(j-1)+2)];
                    end
                end
                
                reps(i) = length(mean_rates);
                
                if reps(i) > 1
%                     odds = [odds; mean_rates(1)];
%                     evens = [evens; mean_rates(2)];

                    odds = [odds; mean(stim_odds)];
                    evens = [evens; mean(stim_evens)];
                end
                
            end
            
            if max(reps)==1
                error('No stimulus repetitions');
            end
            
            explainable_r = corr(odds,evens,'rows','complete');
            
            if nargout < 1
                figure();
                subplot(1,2,1);
                scatter(odds,evens,'jitter','on','jitterAmount',2);
                xlim([0 max(max(odds),max(evens))]);
                ylim([0 max(max(odds),max(evens))]);
                axis square;
                title(sprintf('Odd and Even Reps (r^2 = %.2f)',explainable_r^2));

                subplot(1,2,2);
                hist(reps,1:10);
                xlim([0 11]);
                title('Reps Histogram');
                boxtitle([PFUtil.experName(pf) ' Explainable Variance Analysis']);
            end
            
        end
        
        function [] = plotExplainableVarianceOverTime(pf)
        % PLOTEXPLAINABLEVARIANCEOVERTIME shows how explainable variance
        % changes at different latencies and winsizes for a given p2m file
        % PF
            offsets = 0:10:200;
            winsizes = [30,50,100];
            rs = zeros(length(winsizes),length(offsets));
            for i=1:length(winsizes)
                for j=1:length(offsets)
                    rs(i,j) = RasterUtil.explainableVariance(pf,offsets(j),winsizes(i));
                end
            end
            
            figure();
            plot(offsets,rs','-o');
            title([PFUtil.experName(pf) ' Explainable Variance']);
            ylabel('explainable variance (r^2)');
            xlabel('offset time (ms)');
            ylim([0 1]);
            legend(cellstr(num2str(winsizes','%d ms winsize')));
        end
        
        function [rates_mean,times,rates_std] = psth(pf,pre,post)
        % PSTH computes post-stimulus time histogram (like phist)
        % It returns three vectors: RATES_MEAN, TIMES, and RATES_STD. These
        % If no output arguments are detected then it plots these in a
        % format identical to the one used in phist.
            if nargin < 3
                pre = -200;
                post = 400;
            end
            
            [rates,times] = RasterUtil.rates(pf,pre,post);
            rates_mean = nanmean(rates,1);
            rates_std = nanstd(rates,1);
            
            if nargout < 1
                eshade(times,rates_mean,rates_std);
                hold on;
                plot(times,rates_mean,'-r');
                y_max = max(rates_mean)*1.2;
                line([0 0],[0 y_max],'Color','k');
                ylim([0 y_max]);
                ylabel('s/s');
                xlabel('time (ms)');
            end
        end
        
        function lat = latency(pf)
        % LATENCY estimates a cell's latency LAT (in ms)
        % finds first peak, then chooses the point that lies in between the
        % maxima of the first and second derivatives of the psth. the exact
        % location of this "in-between" point is a weighted average of the
        % maxima (more weight given to the second derivative)
            
            [rates_matrix,ts] = RasterUtil.rates(pf);
            stim_freq = PFUtil.stimulusCarrierFrequency(pf);
            stim_period = (1/stim_freq)*1000;
            
            % average the rates into psth
            rates = mean(rates_matrix,1);
            
            % only examine response from 50-150ms after onset
            MIN_T = 40;
            MAX_T = 150;
            ts_min_i = find(ts==MIN_T,1);
            ts_max_i = find(ts==MAX_T);
            rates_slice = rates(ts_min_i:ts_max_i);
            ts_slice = ts(ts_min_i:ts_max_i);
            
            [peaks,locs] = findpeaks(rates_slice);%,'MINPEAKDISTANCE',round(0.75*stim_period));
            
            % find minimum in first part of rates slice, and make sure that 
            % peak occurs after minimum
            [~,min_i] = min(rates_slice(1:50));
            slice_minimum_t = ts_slice(min_i);
            if length(locs) > 1
                locs = locs(ts_slice(locs) > slice_minimum_t);
            end
            
            % take smallest local maximum that occurs after the minimum as
            % the transient onset peak
            max_i = min(locs);
            max_t = ts_slice(max_i);
            transient_peak_i = ts_min_i + max_i - 1; % adjust back to indices of rates vector
            
            % latency = max_t - 10;
            
            % examine first and second derivatives in the PRE_PEAK_WINDOW
            % milliseconds before the peak
            PRE_PEAK_WINDOW = 40;
            if ts(transient_peak_i - PRE_PEAK_WINDOW) > slice_minimum_t
                PRE_PEAK_WINDOW = ts(transient_peak_i) - slice_minimum_t;
            end
            null_rates = rates((transient_peak_i - PRE_PEAK_WINDOW):transient_peak_i);
            null_ts = ts((transient_peak_i - PRE_PEAK_WINDOW):transient_peak_i);
            null_diff = diff(null_rates);
            null_diff2 = diff(null_diff);
            
            [~,null_diff_peak] = max(null_diff);
            [~,null_diff2_peak] = max(null_diff2);
            
            % take weighted average between two derivative maxima
            lat = null_ts(ceil((null_diff_peak + 2*null_diff2_peak)/3));
        end
        
         function [] = plotRateClusters(pf, n_clusters)
             
            if nargin < 2
                n_clusters = 'auto';
            end
             
            [rates_full,ts_full,triggers] = RasterUtil.rates(pf);
            
            latency = RasterUtil.latency(pf);
            PRE_LATENCY_WINDOW = 20;
            T_MIN = latency - PRE_LATENCY_WINDOW;
            T_MAX_DEFAULT = 250;
            
            stim_period = (1 / PFUtil.stimulusCarrierFrequency(pf))*1000;
            
            T_MAX_COMPUTED = round(T_MIN + stim_period + PRE_LATENCY_WINDOW);
            
            T_MAX = min(T_MAX_DEFAULT,T_MAX_COMPUTED);
            
            rate_i_min = find(ts_full==T_MIN);
            rate_i_max = find(ts_full==T_MAX);
            
            rates = rates_full(:,rate_i_min:rate_i_max);
            rates(isnan(rates)) = 0;
            ts = ts_full(rate_i_min:rate_i_max);
            
            mean_rate = mean(rates(:));
            for i=1:size(rates,1)
                rates(i,:) = rates(i,:) - mean(rates(i,:));% - mean_rate;
                rates(i,:) = rates(i,:) / std(rates(i,:));
            end
            
            rates = rates(~isnan(mean(rates,2)),:);
% 
%             path = ['/lab/results/batch_pstrfs/' stimulus '/'];
%             files_query = sprintf('%s*%.2f*.mat',path,lambda);
%             files = dir(files_query);
% 
%             if strcmp(stimulus,'angleplay')
%                 N_LAGS = 24;
%             else
%                 N_LAGS = 25;
%             end
% 
%             vectors = zeros(length(files),N_LAGS);
%             lags = [];
%             for i=1:length(files)
%                 filename = [path files(i).name];
%                 load(filename);
% 
%                 P = pstrf_result;
%                 lags = P.lags(1:N_LAGS);
%                 vector_sig_acc = zeros(N_LAGS,1);
%                 n_sigs = 0;
%                 for j=1:length(P.eigenvalues_significance)
%                     if P.eigenvalues_significance(j,1) == 1
%                        vector_sig_acc = vector_sig_acc + (P.eigen_vectors(:,j)*P.eigenvalues(j));
%                        n_sigs = n_sigs + 1;
%                     end
%                 end
% 
%                 weighted_eigen_vector = vector_sig_acc / n_sigs;
%                 first_eigen_vector = P.eigen_vectors(:,1);
%                 mean_beta = P.theta_mean;
% 
%                 % Default Scaling (i.e. NONE!)
%                 vectors(i,:) = mean_beta;
% 
%                 % Scale Max to 1
%                 % vectors(i,:) = vectors(i,:) / max(abs(vectors(i,:)));
% 
%                 % Scale to Unit Variance
%                 vectors(i,:) = vectors(i,:) / std(vectors(i,:));
%             end
% 
%             % get rid of nan vectors
%             vectors = vectors(sum(isnan(vectors),2)==0,:);

            % K-means clustering
            if strcmp(n_clusters,'auto')==1
                ks = 2:6;
                silhouette_means = zeros(length(ks),1);
                for i=1:length(ks)
                    cluster_ix = kmeans(rates, ks(i),'replicates',5);
                    silh = silhouette(rates,cluster_ix);
                    silhouette_means(i) = mean(silh);
                end
                [best_silh,best_k_i] = max(silhouette_means);
                n_clusters = ks(best_k_i);
                fprintf('best k=%d. mean silhouette = %.2f\n',n_clusters,best_silh);
            end
            
            [cluster_ix,clusters] = kmeans(rates, n_clusters,'replicates',5);
            
            N_SUBPLOTS = n_clusters + 4;
            N_COLS = 2;
            N_ROWS = ceil(N_SUBPLOTS/N_COLS);
            
            color_palette = cell(6,1);
            color_palette{1} = 'r';
            color_palette{2} = 'b';
            color_palette{3} = 'g';
            color_palette{4} = 'c';
            color_palette{5} = 'm';
            color_palette{6} = [0.5 0.5 0.5];

            % MDS
            dissimilarities = pdist(rates);
            
            % sstress criterion prevents co-location
            [MDS_Y,stress] = mdscale(dissimilarities,2,'criterion','metricsstress');

            figure();
            
            rates_full_mean = nanmean(rates_full,1);
            max_cluster_y = max(clusters(:))*1.2;
            min_psth_y = min(rates_full_mean);
            max_psth_y = 1.2*max(rates_full_mean);
            
            subplot(N_ROWS,N_COLS,1);
           
            rectangle('Position',[min(ts),0,(max(ts)-min(ts)),max_psth_y],...
                'FaceColor',[0.93 0.97 1],...
                'LineStyle','none');
            hold on;
            plot(ts_full,rates_full_mean,'Color',[0.3 0.3 0.3],'LineWidth',1.5);
           % hold on;
            latency = RasterUtil.latency(pf);
            plot(ts,rates_full_mean(rate_i_min:rate_i_max),'-k','LineWidth',1.5);
            stim_freq = PFUtil.stimulusCarrierFrequency(pf);
            stim_period = round((1/stim_freq)*1000);
            for i=1:5
                line([(i-1)*stim_period (i-1)*stim_period],[0 max_cluster_y],'Color',[0.6 0.6 0.6]);
            end
            
            lat_i = find(ts_full==latency);
            mean_rates = mean(rates_full,1);
            marker_y_offset = max_psth_y * 0.05;
            scatter([latency],[mean_rates(lat_i)-marker_y_offset],60,'^','MarkerFaceColor','r','MarkerEdgeColor','none');
            %line([latency latency],[0 max_cluster_y],'Color',[ 0.95 0.95 0.95],'LineStyle','-');
            ylim([min_psth_y max_psth_y]);
            xlim([-100 300]);
            set(gca,'Color','w');%[0.9 0.9 0.9]);
            title('Full PSTH');
            
            subplot(N_ROWS,N_COLS,2);

            for i=1:n_clusters
                scatter(MDS_Y((cluster_ix==i),1),MDS_Y((cluster_ix==i),2),36,color_palette{i});
                hold on;
            end
            title('Multidimensional Scaling');
            
            % Mean Plot
            subplot(N_ROWS,N_COLS,3);
            line([min(ts) max(ts)],[0 0],'Color','k');
            hold on;
            
            % bin ts and rates vectors so that there aren't as many error
            % bars
            ts_binned = SuperUtil.binVector(ts,10);
            rates_mean_binned = SuperUtil.binVector(mean(rates,1),10);
            rates_std_binned = SuperUtil.binVector(std(rates,1),10);
            errorbar(ts_binned,rates_mean_binned,rates_std_binned,'Color',[0.8 0.8 0.8]);
            plot(ts,mean(rates,1),'LineWidth',1.5,'Color','k');
            xlim([min(ts) max(ts)]);
            title('Mean PSTH');
            xlabel('lag (ms)');
            ylabel('rate (s/s)');
            ylim([-max_cluster_y max_cluster_y]);
            
            subplot(N_ROWS,N_COLS,4);
            silhouette(rates,cluster_ix);
            title('Cluster Silhouettes');
            
            
            for i=1:n_clusters
                subplot(N_ROWS,N_COLS,i+4);
                line([min(ts) max(ts)],[0 0],'Color','k');
                hold on;
                plot(ts,clusters(i,:),'LineWidth',1.5,'Color',color_palette{i});
                xlim([min(ts) max(ts)]);
                ylim([-max_cluster_y max_cluster_y]);%ylim([-max_cluster_y max_cluster_y]);
                xlabel('lag (ms)');
                ylabel('rate (s/s)');
                title(sprintf('K-Means Centroid #%d (n=%d)',i,sum(cluster_ix==i)));
            end
            
            if ~isempty(best_silh)
                boxtitle(sprintf('%s Post-Stimulus Response Clusters (silhouette=%.2f)',PFUtil.experName(pf),best_silh));
            else
                boxtitle(sprintf('%s Post-Stimulus Response Clusters',PFUtil.experName(pf)));
            end
         end
        
         function [] = plotPrefNonprefPSTH(pf,transient_min_t,transient_max_t)
         % PLOTPREFNONPREFPSTH compares psths of preferred and
         % non-preferred stimuli
         %
         % The responses to each unique stimulus are sorted based on the
         % intensity of the transient response, defined as the mean rate
         % between TRANSIENT_MIN_T and TRANSIENT_MAX_T.
         %
         % The plots are centered around the null/spontaneous firing rate,
         % and normalized to have a max of 1.

            PLOT_XLIM = [0 200];

            if nargin < 3
                transient_min_t = 60;
                transient_max_t = 100;
            end

            [rates,times] = RasterUtil.rates(pf);
            
            % Find transient firing rates (used for ranking responses)
            transient_min_i = find(times==transient_min_t);
            transient_max_i = find(times==transient_max_t);
            transient_rates = rates(:,transient_min_i:transient_max_i);
            transient_rates_mean = mean(transient_rates,2);
            
            % Threshold set at z-score >= 1
            thresh = mean(transient_rates_mean)+std(transient_rates_mean);
            
            % Preferred rates
            pref_i = transient_rates_mean >= thresh;
            pref_rates = rates(pref_i,:);
            n_pref = size(pref_rates,1);
            
            % Non-preferred rates
            nonpref_i = ~pref_i;
            nonpref_rates = rates(nonpref_i,:);
            n_nonpref = size(nonpref_rates,1);
            
            % Compute PSTHs
            psth_mean = mean(rates,1);
            psth_pref = mean(pref_rates,1);
            psth_nonpref = mean(nonpref_rates,1);
                        
            % Slice PSTH to only include data in PLOT_XLIM bounds
            plot_min_i = find(times==min(PLOT_XLIM));
            plot_max_i = find(times==max(PLOT_XLIM));
            times_slice = times(plot_min_i:plot_max_i);
            psth_mean = psth_mean(plot_min_i:plot_max_i);
            psth_pref = psth_pref(plot_min_i:plot_max_i);
            psth_nonpref = psth_nonpref(plot_min_i:plot_max_i);
            
            % Find null firing rate (in first 40ms)
            null_rate_mean = mean(psth_mean(1:45));
            
            % Normalize PSTHs
            psth_mean_n = psth_mean - null_rate_mean;
            psth_mean_n = psth_mean_n / max(psth_mean_n);
            psth_pref_n = psth_pref - null_rate_mean;
            psth_pref_n = psth_pref_n / max(psth_pref_n);
            psth_nonpref_n = psth_nonpref - null_rate_mean;
            psth_nonpref_n = psth_nonpref_n / max(psth_nonpref_n);
            
            % Generate Figure
            figure();
            subplot(2,1,1);
            plot(times_slice,psth_mean,'-k');
            hold on;
            plot(times_slice,psth_pref,'-b');
            plot(times_slice,psth_nonpref,'-r');
            xlim(PLOT_XLIM);
            legend('mean',['preferred (n=' num2str(n_pref) ')'],['non-preferred (n=' num2str(n_nonpref) ')']);
            title([PFUtil.experName(pf) ' Temporal Dynamics']);
            ylabel('s/s');
            xlabel('time (ms)');
            
            subplot(2,1,2);
            line(PLOT_XLIM,[0 0],'Color','k');%plot(times_slice,psth_mean_n,'-k');
            hold on;
            p1 = plot(times_slice,psth_pref_n,'-b');
            p2 = plot(times_slice,psth_nonpref_n,'-r');
            xlim(PLOT_XLIM);
            ylim([-1.1 1.1]);
            legend([p1,p2],['preferred (n=' num2str(n_pref) ')'],['non-preferred (n=' num2str(n_nonpref) ')']);
            title('Normalized Comparison');
            ylabel('normalized response');
            xlabel('time (ms)');
            box on;
         end
         
         function [] = plotExperPSTH(pf_list)
         % PLOTEXPERPSTH plots multiple psths from the same experiment
         % together
         %
         % PF_LIST is a cell array with multiple p2m files
         

            PLOT_XLIM = [0 250];

            times = zeros(range(PLOT_XLIM)-1,1);
            psths = zeros(length(pf_list),range(PLOT_XLIM)-1);
            psths_std = zeros(length(pf_list),range(PLOT_XLIM)-1);
            psths_n = zeros(length(pf_list),range(PLOT_XLIM)-1);
            psths_std_n = zeros(length(pf_list),range(PLOT_XLIM)-1);
            task_names = cell(length(pf_list),1);
            
            % Get PSTHs for each p2m file in list
            for i=1:length(pf_list)
                [psths(i,:),times,psths_std(i,:)] = RasterUtil.psth(pf_list{i},min(PLOT_XLIM),max(PLOT_XLIM));
                task_names{i} = PFUtil.taskName(pf_list{i});
            end
            
            % Calculate spontaneous / null rate (in first 45ms)
            null_rates = mean(psths(:,1:45),2);
                                    
            
            % Normalize PSTHs
            for i=1:length(pf_list)
                psths_n(i,:) = psths(i,:) - null_rates(i);
                max_rate = max(psths_n(i,:));
                psths_n(i,:) = psths_n(i,:) / max_rate;
                psths_std_n(i,:) = psths_std(i,:) / max_rate;
            end
            
            % Generate Figure
            figure();
            subplot(2,1,1);
            plot(times,psths');
            xlim(PLOT_XLIM);
            title([PFUtil.experName(pf_list{1}) ' Temporal Dynamics']);
            ylabel('s/s');
            xlabel('time (ms)');
            legend(task_names);
            
            subplot(2,1,2);
            plot(times,psths_n');
            xlim(PLOT_XLIM);
            ylim([-1.1 1.1]);
            title('Normalized Comparison');
            ylabel('normalized response');
            xlabel('time (ms)');
            box on;
            legend(task_names);
            line(PLOT_XLIM,[0 0],'Color','k');
         end
    end
end