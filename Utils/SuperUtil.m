classdef SuperUtil < handle
    % SUPERUTIL
    
    properties
             
    end
    
    properties (Constant)
        
        DT = 0.001;
        
    end
    
    methods

        function SU = SuperUtil()
        end
        
    end
    
    methods(Static)
        
        function [ngrams, ngrams_error] = ngramTransitions(pf, max_isi, bootstrap)
        % NGRAMTRANSITIONS
        %
        % Params:
        % - pf
        % - max_isi
        % - bootstrap: (0/1) if ==1, call recursively to do bootstrap
            
            N_MAX_PRIORS = 10;
            
            if nargin < 3
                bootstrap = 0;
            end
            
            if bootstrap==1
                N_DRAWS = 100;
                N_RECS = length(pf.rec);
                ngrams_matrix = zeros([N_DRAWS N_MAX_PRIORS+2]);
                
                for i=1:N_DRAWS
                    new_pf = pf;
                    new_pf.rec = new_pf.rec(randi(N_RECS,N_RECS,1));
                    res = SuperUtil.ngramTransitions(new_pf,max_isi);
                    ngrams_matrix(i,:) = res;
                end
                ngrams = mean(ngrams_matrix,1);
                ngrams_error = std(ngrams_matrix,1);
                return;
            end
                    
            priors = [];
            posts = [];
                        
            for n=1:length(pf.rec)
                spikes = pf.rec(n).spike_times;
                
                % this is a shortcut! there must be a better way to decide
                % whether to count this data or not... like looking for
                % stimulate_start and stimulate_end?
                if isempty(spikes)
                    continue;
                end
                
                % Find the number of prior spikes and if there was a spike
                % afterward for each spike (i.e. the priors and posts
                % vectors)
                [rec_priors, rec_posts] = SuperUtil.spikes2priors(spikes, max_isi);
                priors = [priors; rec_priors];
                posts = [posts; rec_posts];
            end
                        
            ngrams = zeros([1 N_MAX_PRIORS+1]);
            
            % 0th probability
            ngrams(1) = length(priors) / SuperUtil.duration(pf);
            
            % 1st...nth order:
            for n=0:N_MAX_PRIORS
                ngrams(n+2) = sum((priors >= n) & (posts==1)) / sum(priors >= n);
            end
            
            % replace nans with 0s
            ngrams(isnan(ngrams)) = 0;
            
            % Observation: my intuition was that as max_isi increases, the
            % probabilities will necessarily get larger. However this is
            % not always the case: although the total number of spikes with
            % priors will always increase, the number of spikes-with-priors
            % that are followed by another spike may not.
        end
        
        function [priors,posts] = spikes2priors(spikes, max_isi)
        % SPIKES2PRIORS takes a standard p2m spike train and returns two
        % vectors (of the same length): PRIORS is the number of spikes that
        % precede each spike (where each of these spikes is separated by
        % MAX_ISI), and POSTS is a binary vector reporting whether another
        % spike within MAX_ISI follows each spike.
            
            priors = zeros(length(spikes), 1);
            
            isis = diff(spikes);
            isis(length(spikes))=100000; % hack: added to prevent segfault 
                                         % in step 2 below, since diff
                                         % reduces vector length by 1
                                         
            % Step 1: count how many spikes preceeded each spike
            % (aka the priors)

            % find indices of bursters
            bursters = isis <= max_isi;
            bursters_is = find(bursters==1) + 1; % add 1 since we're adding prior count to next spike
            
            % loop over bursters
            for i=1:length(bursters_is)
                % if the last spike was within MAX_ISI ms, look at the 
                % priors count of the previous spike and add 1 to it 
                % #pseudorecursion
                priors(bursters_is(i)) = (priors(bursters_is(i)-1) + 1);
            end
            
            % Step 2: note whether another spike follows this one
            posts = isis <= max_isi;
        end

        function [] = plotNgramPoisson(pf)
        % PLOTNGRAMPOISSON 
            N_DRAWS = 2;
            WINSIZES = [10,20,50];
            MAX_ISI = 4;
            [ngrams,ngrams_error] = SuperUtil.ngramTransitions(pf, MAX_ISI,1);

            APW = AnglePlayWriter(pf);

            ngrams_p = zeros(N_DRAWS,length(WINSIZES),length(ngrams));

            for i=1:N_DRAWS
                for j=1:length(WINSIZES)
                    fprintf('winsize=%d, draw=%d\n',WINSIZES(j),i);
                    new_pf = APW.poissonSpikes(WINSIZES(j));
                    new_pf.src = [new_pf.src '-' num2str(WINSIZES(j))];
                    ngrams_p(i,j,:) = SuperUtil.ngramTransitions(new_pf, MAX_ISI); % only plots on first iteration
                end
            end

            % calculate average
            ngrams_p_mean = squeeze(mean(ngrams_p,1));
            
            % Plot
            h = errorbar(1:12,ngrams,ngrams_error,'-k','LineWidth',2);

            hold on;
            plot(ngrams_p_mean','-o');
            legend('observed','10ms bins','20ms bins','50ms bins');
            set(gca,'YGrid','on');
            ylim([0 1]);
            xlim([1 11]);
            set(gca,'XTick',1:12,'XTickLabel',{'p(1)','p(2|1)','p(3|2)','p(4|3)', 'p(5|4)','p(6|5)','p(7|6)','p(8|7)','p(9|8)','p(10|9)','p(11|10)'});

        end
        
        function [] = plotNgrams(pf,winsize)
        % PLOTNGRAMS plots the combined spike transition probability data
        % visualization. This includes a PSTH, a spikes per bin histogram,
        % and both the observed and poisson ngram probabilties.
        %
        % Params:
        % - pf: p2m file
        % - winsize: sliding window size used for calculating histogram
            
            [counts,n] = SuperUtil.spikesPerBinHist(pf,winsize);
            counts_n = counts / sum(counts);
            ngrams_3 = SuperUtil.ngramTransitions(pf,3);
            ngrams_4 = SuperUtil.ngramTransitions(pf,4);
            ngrams_5 = SuperUtil.ngramTransitions(pf,5);
            
            
            figure();

            plot_x_range = [0 10];

            subplot(3,2,1);
            psth = phist(pf);
            eshade(psth(:,1),psth(:,2),psth(:,3));
            hold on;
            plot(psth(:,1),psth(:,2),'-r');
            title('PSTH');
            ylim([0 (max(psth(:,2)) + max(psth(:,3)))]);


            subplot(3,2,2);
            bar(n, counts_n);%, '-o');
            title(['Spikes per Bin Histogram (bin = ' num2str(winsize) 'ms)']);
            xlabel('n_{spikes}');
            ylabel('fraction of bins');
            xlim(plot_x_range - 0.5);
            ylim([0 max(counts_n)*1.05]);


            subplot(3,2,3:4);
            plot(ngrams_3,'-or');
            hold on;
            plot(ngrams_4,'-ok','LineWidth',2);
            plot(ngrams_5,'-ob');
            title('Transition Probabilities');
            set(gca,'XTick',1:11,'XTickLabel',{'p(1)','p(2|1)','p(3|2)','p(4|3)', 'p(5|4)','p(6|5)','p(7|6)','p(8|7)','p(9|8)','p(10|9)','p(11|10)'});
            ylim([0 1]);
            xlim([1 11]);
            set(gca,'YGrid','on');
            legend('ISI_{max} = 3ms','ISI_{max} = 4ms','ISI_{max} = 5ms');
            
            subplot(3,2,5:6);
            SuperUtil.plotNgramPoisson(pf);
            title('Observed vs. Inhomogeneous Poisson-Generated Transition Probabilities (ISI_{max}=4ms)');
            
            boxtitle([pf.src ' winsize=' num2str(winsize)]);
        end
        
        function [count,centers] = plotISI(pf)
        % PLOTISI plots the inter-spike interval distribution
            isis = [];
            
            xmax = 200;
            
            for i=1:length(pf.rec)
                spikes = pf.rec(i).spike_times;
                dts = diff(spikes);
                isis = [isis; dts];
            end
            
            figure();
            [count,centers] = hist(isis, 0:1:xmax);
            bar(centers,(count / sum(count(:))));
            xlim([0 100]);
            
        end
        
        function [isi_matrix_n] = burstISI(pf, MAX_BURST)
            
            BURST_SIZES = 1:MAX_BURST;
            MAX_BURST_ISI = 5;
            MAX_ISI = 200;
            BUFFER_SIZE = length(pf.rec) * 500;
            
            isi_matrix = zeros(length(BURST_SIZES), MAX_ISI);
            isis_raw = zeros(length(BURST_SIZES), BUFFER_SIZE);
            isis_raw_max_i = zeros(length(BURST_SIZES),1);
            
            for r=1:length(pf.rec)
                
                spikes = pf.rec(r).spike_times;
                
                if isempty(spikes)
                    continue;
                end
                
                isis = diff(spikes);
                [priors,~] = SuperUtil.spikes2priors(spikes,MAX_BURST_ISI);

                for i=1:length(BURST_SIZES)
                    
                    bursters_is = priors >= (BURST_SIZES(i)-1);
                    % spikes_bursting = spikes(bursters_is);
                    isis_bursting = isis(bursters_is(1:length(isis)));
                    start_i = isis_raw_max_i(i) + 1;
                    stop_i = start_i + length(isis_bursting) - 1;
                    isis_raw(i,start_i:stop_i) = isis_bursting;
                    isis_raw_max_i(i) = stop_i;
                end
            end
            
            for i=1:length(BURST_SIZES)
                isi_matrix(i,:) = hist(isis_raw(i,1:isis_raw_max_i(i)),1:MAX_ISI);
            end
            
            % Normalize ISI Matrix
            sums = sum(isi_matrix,2);
            sums_matrix = repmat(sums,[1 MAX_ISI]);
            isi_matrix_n = isi_matrix ./ sums_matrix;
        end
        
        function [] = plotBurstISI(pf, MAX_BURST)
            
            BURST_SIZES = 1:MAX_BURST;
            MAX_ISI_TO_DISPLAY = 50;
            
            isi_matrix_n = SuperUtil.burstISI(pf, MAX_BURST);
            
            figure();
            for i=1:length(BURST_SIZES)
                subplot(length(BURST_SIZES),1,BURST_SIZES(i));
                bar(1:length(isi_matrix_n(i,:)),isi_matrix_n(i,:));
                xlim([0 MAX_ISI_TO_DISPLAY]);
                ylim([0 max(isi_matrix_n(:))]);
                title([num2str(i-1) ' Prior Spikes']);

                if i==length(BURST_SIZES)
                    xlabel('ISI (ms)');
                end
            end
            boxtitle([SuperUtil.experName(pf) ' Burst ISI']);
        end
        
        function [x,ys,errs] = plotBootstrappedBurstISI(pf, MAX_BURST)
            
            N_DRAWS = 5;
            BURST_SIZES = 1:MAX_BURST;
            MAX_ISI_TO_DISPLAY = 30;
            
            [isis_mean,isis_std] = PFUtil.bootstrap(pf,@SuperUtil.burstISI, N_DRAWS, MAX_BURST);
            
            x = 1:MAX_ISI_TO_DISPLAY;
            ys = isis_mean(:,1:MAX_ISI_TO_DISPLAY);
            errs = isis_std(:,1:MAX_ISI_TO_DISPLAY);
            
            figure();
            for i=1:length(BURST_SIZES)
                subplot(length(BURST_SIZES),2,2*(BURST_SIZES(i)-1)+1);
                
                errorbar(x,ys(i,:),errs(i,:));
                xlim([0 MAX_ISI_TO_DISPLAY]);
                ylim([0 max(ys(:))*1.1]);
                title([num2str(i-1) ' Prior Spikes']);

                if i==length(BURST_SIZES)
                    xlabel('ISI (ms)');
                end
            end
            boxtitle([SuperUtil.experName(pf) ' Burst ISI']);
            
            subplot(length(BURST_SIZES),2,2:2:(length(BURST_SIZES)*2));
            cdfs = cumsum(ys,2);
            plot(x,cdfs);
            legendCell = cellstr(num2str((BURST_SIZES-1)', '%d Prior Spikes'));
            legend(legendCell,'Location','SouthEast');
            title('ISI CDFs');
            grid on;
            xlim([0 MAX_ISI_TO_DISPLAY]);
            ylim([0 1]);
            
        end
        
        function [counts,rates] = rateHist(pf, winsize)
            
            if nargin < 2
                winsize = 20;
            end
            
            rates_acc = [];
            
            for i=1:length(pf.rec)
                
                spikes = pf.rec(i).spike_times;
                [~,r] = SuperUtil.spikes2rates(spikes,'winsize',winsize);
                rates_acc = [rates_acc; r];
                
            end
            
            [raw_counts,raw_rates] = hist(rates_acc,20);
            
            % get rid of points without samples. the nature of the spike
            % histogram procedure means that some bins will be zero.
            counts = raw_counts(raw_counts > 0);
            rates = raw_rates(raw_counts > 0);
        end
        
        function n = spikeCount(pf)
            warning('Deprecated. Use PFUtil.spikeCount() instead');
            n = PFUtil.spikeCount(pf);
        end
        
        function d = duration(pf)
        % DURATION finds an approximation of the total time of the
        % experiment across all trials. It is not precise, because it finds
        % the start/stop time by looking for the times of the first and
        % last spikes.
            d = 0;
            for i=1:length(pf.rec)
                spikes = pf.rec(i).spike_times;
                d = d + (max(spikes) - min(spikes));
            end
        end
        
        function [hist_n, hist_bins] = autocorrelogram(pf, BIN_WIDTH)
            
            SEARCH_WINDOW = 2000; % measure correlations with +/- WINDOW elements in spikes vector
            PLOT_WINDOW = 250; % xlim = +/- PLOT_WINDOW
            
            if nargin < 2
                BIN_WIDTH = 3;
            end
            
            
            spikes = PFUtil.concatenateSpikes(pf);
            n_spikes = length(spikes);
            
            intervals = nan((2*SEARCH_WINDOW + 1) * n_spikes, 1);
            max_int_i = 0;
            
            for i=1:length(spikes)
                min_i = max(1, i - SEARCH_WINDOW);
                max_i = min(n_spikes, i + SEARCH_WINDOW);
                spike_intervals = spikes(min_i:max_i) - spikes(i);
                
                intervals_is = (max_int_i+1):(max_int_i+length(spike_intervals));
                intervals(intervals_is) = spike_intervals;
                max_int_i = max(intervals_is);
            end
            
            intervals = intervals(~isnan(intervals));
            
            % filter out intervals outside of PLOT_WINDOW
            intervals = intervals(abs(intervals) <= PLOT_WINDOW);
            
            
            hist_bins = -PLOT_WINDOW:BIN_WIDTH:PLOT_WINDOW;
            hist_n = hist(intervals,hist_bins);
            
            % remove center
            hist_n(floor(length(hist_n)/2)+1)=0;

            if nargout < 1
                figure();
                bar(hist_bins,hist_n);
                xlim([-PLOT_WINDOW, PLOT_WINDOW]);
                title(sprintf('%s  Autocorrelogram (%d ms bins)',PFUtil.experName(pf),BIN_WIDTH));
            end
        end
        
        function [hist_n, hist_bins] = bootstrappedAutocorrelogram(pf)
            
            N_DRAWS = 50;
            N_RECS = length(pf.rec);
            [hist_n_test,hist_bins] = SuperUtil.autocorrelogram(pf);
            
            hist_n_acc = zeros(size(hist_n_test));
            
            for i=1:N_DRAWS
                new_pf = pf;
                new_pf.rec = new_pf.rec(randi(N_RECS,N_RECS));
                [draw_n,~] = SuperUtil.autocorrelogram(new_pf);
                hist_n_acc = hist_n_acc + draw_n;
            end
            
            hist_n = hist_n_acc / N_DRAWS;
        end
        
        function t = transientOnsetTime(pf)
            
            % Plan:
            % 1) get psth from phist()
            % 2) Measure mean and std of firing rates between -50 and 50ms
            % 3) Find max firing rate between 50 and 150
            % 4) Go backwards and find the first point that is greater than
            %    the null mean + std.
            
            psth = phist(pf);
            
            % Find null mean and std
            null_t_min = find(psth(:,1)==0);
            null_t_max = find(psth(:,1)==40);
            null_mean = mean(psth(null_t_min:null_t_max,2));
            null_std = std(psth(null_t_min:null_t_max,2));
            
            % Find transient peak
            transient_t_min = find(psth(:,1)==50);
            transient_t_max = find(psth(:,1)==150);
            [peak_fr, transient_peak] = max(psth(transient_t_min:transient_t_max,2));
            transient_peak = transient_peak + transient_t_min - 1;
            
            fr = peak_fr;
            i = transient_peak;
            while fr >= null_mean
                fr = psth(i,2);
                i = i - 1;
            end
            
            t = psth(i,1);
            
            phist(pf);
            hold on;
            psth_ylim = ylim;
            line([t t],psth_ylim,'Color','b');
            xlim([-100 200]);
            
        end
        
        function r = meanRate(pf)
        % MEANRATE calculates the mean rate from all trials
            
            T_START = -100;
            T_STOP = 100;
            
            rast = prast(pf);
            
            t_start_i = find(rast.time==T_START);
            t_stop_i = find(rast.time==T_STOP);
            
            spikes = rast.data(:,t_start_i:t_stop_i);
            
            r = nanmean(spikes(:)) * 1000;
        end
                
        function [ns,counts,bin] = poissonBinHist(pf, meanRate)
            
            INTERVAL = 0.5; % also try INTERVAL = 2, N_MAX = 20
            N_MIN = 0;
            N_MAX = 5;
            
            if nargin < 2
                meanRate = SuperUtil.meanRate(pf);
            end
            
            bin = (1 / meanRate)*1000;
            ns = N_MIN:INTERVAL:N_MAX;
            winsizes = zeros(size(ns));
            counts = zeros(size(ns));
            means = zeros(size(ns));
            
            parfor i=1:length(ns)
                winsizes(i) = round(ns(i) * bin);
                [hist_counts,hist_n] = SuperUtil.spikesPerBinHist(pf,winsizes(i));
                counts(i) = sum(hist_counts(:));%(hist_n==3 | hist_n==4 | hist_n==5));
                means(i) = sum(hist_counts .* hist_n) / sum(hist_counts);
            end
            
            counts = means;
        end
        
        function [] = plotPoissonBinHist(pf, winsize)
            
            N_DRAWS = 12;
            if nargin < 2
                POISSON_WINDOW = 20;
            else
                POISSON_WINDOW = winsize;
            end
            REFRACTORY_PERIOD = 3;
            N_REPETITIONS = 2;
            
            mean_rate = SuperUtil.meanRate(pf);
            [ns_obs, counts_obs,bin] = SuperUtil.poissonBinHist(pf, mean_rate);
            
            counts_obs_matrix = zeros(N_DRAWS,length(counts_obs));
            
            for i=1:N_DRAWS
                bootpf = pf;
                fprintf('observed draw %d\n',i);
                bootpf.rec = bootpf.rec(randi(length(pf.rec),length(pf.rec),1));
                [~, counts_obs_matrix(i,:)] = SuperUtil.poissonBinHist(bootpf, mean_rate);
            end
            
            APW = AnglePlayWriter(pf);
            fprintf('homo draw\n');
            randpf_homo = APW.homogeneousPoisson(mean_rate, REFRACTORY_PERIOD,'repetitions',N_REPETITIONS);
            [~, counts_rand_homo] = SuperUtil.poissonBinHist(randpf_homo, mean_rate);
            
            randpf_inhomo = APW.inhomogeneousPoisson(POISSON_WINDOW, REFRACTORY_PERIOD,'repetitions',N_REPETITIONS);
            fprintf('inhomo draw\n');
            [~, counts_rand_inhomo] = SuperUtil.poissonBinHist(randpf_inhomo, mean_rate);
                        
            errorbar(ns_obs,mean(counts_obs_matrix,1),std(counts_obs_matrix,1),'-k','LineWidth',1);
            hold on;
            plot(ns_obs,mean(counts_rand_homo,1),'-ro');
            plot(ns_obs,mean(counts_rand_inhomo,1),'-bo');
            set(gca,'XTick',0:2:max(ns_obs));
            %grid on;
            legend('Observed',['Poisson, FR=' num2str(round(mean_rate))], ['Poisson, FR Varying, winsize=' num2str(POISSON_WINDOW)],'Location','SouthEast');
            ylabel('Expected # Spikes Per Bin');
            xlabel(['n (Bin Width = n*\lambda^{-1} = n*' num2str(round(bin)) 'ms)']);
            title([SuperUtil.experName(pf) ' Expected Spikes Per Bin']);
            xlim([min(ns_obs), max(ns_obs)]);
        end
        
        function [counts,n] = spikesPerSingleBinHist(pf,latency,winsize)
            raster = prast(pf);
            
            HIST_MAX = ceil(winsize / 2);
            n = 0:HIST_MAX;
            n_trials = size(raster.data,1);
            ts = latency:latency+winsize;
            
            % find indices of start and stop times
            min_i = find(raster.time==latency);
            max_i = find(raster.time==latency+winsize);
            
            counts_acc = zeros(n_trials, 1);
            firsts = [];
            for i=1:n_trials
                slice = raster.data(i,min_i:max_i);
                counts_acc(i) = sum(slice);
                first_i = find(slice,1,'first');
                if ~isempty(first_i)
                    firsts = [firsts; raster.time(min_i + first_i - 1)]; % minus 1 fixes 0-based indexing bullshit
                end
            end
            counts = hist(counts_acc,n);
            
            %firsts_hist = hist(firsts,ts);
            
            % plot if no output arguments provided, just like hist() does
            if nargout < 1
                figure();
                
                subplot(3,1,1);
                SuperUtil.plotPsthWindow(pf,latency,winsize);
                
                subplot(3,1,2);
                hist(firsts);
                title('Time of First Spike');
                xlabel('ms');
                
                subplot(3,1,3);
                counts_ratio = (counts / sum(counts)) * 100;
                last_nonzero_i = find(counts_ratio,1,'last');
                barh(n(1:last_nonzero_i),counts_ratio(1:last_nonzero_i));
                set(gca,'YDir','reverse');
                title('Spikes Per Bin Histogram');
                xlabel('%');
                ylabel('Spikes In Bin');
                ylim([-0.5 last_nonzero_i - 0.5]);
                boxtitle(sprintf('%s Spikes in %d-%dms', SuperUtil.experName(pf), latency, latency + winsize));
            end
        end
        
        function [counts,n] = spikesPerBinHist(pf, winsize)
        % SPIKESPERBINHIST moves a sliding bin (WINSIZE ms wide) ms-by-ms
        % across each trial in the experiment PF, and creates a histogram of
        % the results. It is useful for measuring the sparsity of a cell's
        % response, and for comparing spike statistics to Poisson
        % predictions.
            
            HIST_MAX = 200;
            n = 0:HIST_MAX;
            
            if winsize==0
                counts = zeros(size(n));
                return;
            end
            
            counts_acc = [];
            
            for i=1:length(pf.rec)
                
                spikes = pf.rec(i).spike_times;
                [~,~,trial_counts] = SuperUtil.spikes2rates(spikes,'winsize',winsize);
                counts_acc = [counts_acc; trial_counts];
                
            end
            
            counts = hist(counts_acc,n);
            
            % plot if no output arguments provided, just like hist() does
            if nargout < 1
                bar(n,counts);
            end
        end
        
        function [ts, rates, counts] = spikes2rates(spikes, varargin)
        % SPIKES2RATES converts a vector of spike times into a vector of
        % firing rates, calculated using a sliding bin of 20ms
        %
        % varagin options:
        %  mint
        %  maxt
        %  winsize
        
            p = inputParser;
            addRequired(p,'spikes');
            addParameter(p,'mint',-1);
            addParameter(p,'maxt',-1);
            addParameter(p,'winsize',20);
            parse(p,spikes,varargin{:});
        
            mint = p.Results.mint;
            maxt = p.Results.maxt;
            winsize = p.Results.winsize;
            dt = winsize / 1000;
            
            if mint == -1
                mint = min(spikes) - round(winsize/2);
                maxt = max(spikes) + round(winsize/2);
            end
            
            ts = mint:1:maxt;
            spike_is = arrayfun(@(x) find(ts == x,1,'first'), spikes );
            
            % convet into binary vector, one element per ms, where a 1 is a
            % spike
            spike_vector = zeros(length(ts),1);
            spike_vector(spike_is) = 1;
            counts = smooth(spike_vector, winsize) * winsize;
            rates = counts / dt;
        end
        
        function [sliced,map] = sliceVector(vector, padding)
            
            % Solution inspired by this stackoverflow answer:
            % http://stackoverflow.com/questions/3274043/finding-islands-of-zeros-in-a-sequence
            vector(vector >= 0.5) = 1;
            vector(vector < 0.5) = 0;
            
            dsig = diff([0 vector 0]);
            startIndex = find(dsig > 0) - padding;
            endIndex = find(dsig < 0)-1 + padding;
            % duration = endIndex-startIndex+1
            
            if startIndex(1) <= 0
                startIndex(1) = 1;
            end
            
            n_islands = length(startIndex);
            map = zeros(n_islands, 4);
            map(:,1) = startIndex;
            map(:,2) = endIndex;
            map(1,3) = 1;
            
            durations = map(:,2) - map(:,1) + ones(n_islands,1);
            
            for i=1:n_islands
                if i == 1
                    map(i,3) = 1;
                else
                    map(i,3) = map(i-1,4) + 1;
                end
                map(i,4) = map(i,3) + durations(i) - 1;
            end
            
            sliced = zeros(sum(durations),1);
            for i=1:n_islands
                sliced(map(i,3):map(i,4)) = vector(map(i,1):map(i,2));
            end
            
        end
        
        function sliced = sliceMatrix(matrix, padding)
            
            n_rows = size(matrix,1);
            sliced = cell(n_rows, 2);
            
            for i=1:n_rows
                [sliced{i,1}, sliced{i,2}] = SuperUtil.sliceVector(matrix(i,:), padding);
            end
        end

        function binned = binVector(vector, bin_width)
            n_bins = floor(length(vector)/bin_width);
            vector_truncated = vector(1:bin_width * n_bins);
            reshaped = reshape(vector_truncated, bin_width, n_bins);
            binned = mean(reshaped,1)';
        end
        
        function binned = binMatrix(matrix, bin_width, dim)
            if dim==2
                n_bins = floor(size(matrix, dim)/bin_width);
                matrix_truncated = matrix(:, 1:n_bins*bin_width);
                matrix_reshaped = reshape(matrix_truncated,[size(matrix,1), bin_width, n_bins]);
                binned = squeeze(mean(matrix_reshaped, 2));
            else
                error('dim=1 not yet implemented');
            end
        end
        
        function [] = plotPsthWindow(pf, latency, winsize)
            phist(pf,'pre',0,'post',300,'suptitle',0);
            hold on;
            y_range = ylim;
            line([latency latency],[min(y_range) max(y_range)],'LineStyle','--','Color','b');
            line([latency+winsize latency+winsize],[min(y_range) max(y_range)],'LineStyle','--','Color','b');
            ylim(y_range);
            title('PSTH','fontweight','bold');
        end
        
        function name = experName(pf)
            warning('Deprecated. Use PFUTil.experName() instead.');
            name = PFUtil.experName(pf);
        end

        
    end
end