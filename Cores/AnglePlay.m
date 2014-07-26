classdef AnglePlay < handle
    %ANGLEPLAY analyzes data from AnglePlay trials
    %   analyzes data collected from curvplay task using angle_* stimuli
    %   set.
    %
    %   @TODO
    %   - orientation tuning curve (collapse across coefficients)
    %   - coefficient tuning curve (collapse across orientation)
    
    properties
        % p2m file
        pf
        
        % curvplay task version number
        task_version
        
        % latency
        latency
        
        % winsize
        winsize
        
        % plot title str
        plot_title
        
        % name of experiment (e.g. romeo0284.curvplay.005)
        exper
        
        % directory where stimulus images are stored
        imagedir
        
        % pype params
        params
        
        % cell matrix w/ blurred stimulus images
        blurredStimulusImages
        
        % length of stimulus image used in experiment (in pixels)
        stimulusSize
        slice_fraction
        
        % stimulus and response data
        responses
        
        max_rate
        
        % null response firing rate stats
        null_mean
        null_std
        
        % stimuli set properties
        N_STIMULI
        orientations
        coefficients
        strokes
        types
        
        
    end
    
    methods
        function AP = AnglePlay(pf, latency, winsize)
        % ANGLEPLAY constructor for AnglePlay class
        %   PARAM pf            (p2m object) data from angleplay trial
        %   PARAM offset        (int) start of bin / cell's latency, in ms
        %   PARAM winsize       (int) size of bin, in ms
        %   RETURN AP           (AnglePlay object)

            AP.pf = getpf(pf);
            AP.latency = latency;
            AP.winsize = winsize;
                
            if isempty(strfind(AP.pf.rec(1).params.X_version_info{2}, 'curvplay'))
                error('not curvplay file');
            end
            
            % name of experiment (e.g. romeo0284.curvplay.005)
            [~,exper_name,exper_no] = fileparts(AP.pf.src);
            AP.exper = [exper_name exper_no];
            AP.plot_title = [AP.exper ' ' num2str(AP.latency) '-' num2str(AP.latency + AP.winsize) 'ms'];
            
            AP.params = AP.pf.rec(1).params;
            AP.imagedir = AP.params.imagedir;
            
            version_str = strjoin(AP.params.X_version_info);
            version_match = regexp(version_str,'curvplay.py\s(\d*)','tokens');
            AP.task_version = str2double(strjoin(version_match{:}));
                        
            AP.responses = AP.findResponses(latency, winsize);
            
            AP.orientations = unique(AP.responses(:,2));
            AP.coefficients = unique(AP.responses(:,4));
            AP.strokes = unique(AP.responses(:,6));
            AP.types = [0, 1]; % parabola=0, abs=1
            
            N_TYPES = length(AP.types); % parabola and abs
            AP.N_STIMULI = length(AP.orientations)*length(AP.coefficients)*length(AP.strokes)*N_TYPES;
            
            % @todo find max rate some other way, findVariance is useless
            v = AP.findVariance();
            AP.max_rate = max(v.means);
            
            [AP.null_mean, AP.null_std] = AP.findAcausalRates();
        end
        
        function [responses] = findResponses(AP, lat, winsize)
            
            % Respdata columns:
            %
            % 1 response        (s/s)
            % 2 ori             (deg)
            % 3 type            (0==parabola, 1==abs)
            % 4 coefficient
            % 5 scale           (mult)
            % 6 stroke width
            % 7 image file #    (number of stim's png file)
            % 8 inum            (stim's index in IMAGE_INFO [zero-based])
            
            opt.lat = lat;
            opt.winsize = winsize;
            opt.unit = 'ttl';			% default to ttl
                        
            image_info = AP.pf.rec(1).params.IMAGE_INFO;
            
            frame_i = 1;
            responses = zeros(100000,8);
            for n = 1:length(AP.pf.rec)
                if AP.pf.rec(n).result(1) == 'C'
                    [~, maxt] = PFUtil.findEvents(AP.pf, n, 'eye_stop');
                elseif strcmp(AP.pf.rec(n).result, 'E FixBreak')
                    [~, maxt] = PFUtil.findEvents(AP.pf, n, 'fix_lost');
                else
                    continue;
                end
                [ix, ts] = PFUtil.findEvents(AP.pf, n, 'frameflip');
                
                for evn = 1:length(ix)
                    if ts(evn)+opt.lat+opt.winsize <= maxt
                        
                        % Get firing rate
                        firing_rate = p2mgetspikes(AP.pf, n, ...
                            ts(evn)+lat, ts(evn)+opt.lat+opt.winsize, ...
                            opt.unit);
                        firing_rate = 1000 * length(firing_rate) ./ opt.winsize;
                        trigger = AP.pf.rec(n).ev_e{ix(evn)};
                        stimParams = AnglePlayUtil.trigger2stimulusParams(AP.pf,trigger);
                        
                        responses(frame_i,1) = firing_rate;
                        responses(frame_i,2) = stimParams.ori;
                        responses(frame_i,3) = stimParams.type;
                        responses(frame_i,4) = stimParams.coefficient;
                        responses(frame_i,5) = stimParams.scale;
                        responses(frame_i,6) = stimParams.stroke;
                        responses(frame_i,7) = stimParams.image_file_n;
                        responses(frame_i,8) = stimParams.inum;
                        frame_i = frame_i + 1;
                    end
                end
            end
            % trim respdata
            responses = responses(1:(frame_i-1),:);
            
        end
        
        function [null_mean,null_std] = findAcausalRates(AP)
            
            % find responses at acausal latencies
            latencies = -200:30:-30;
            
            rates = [];
            for lat=latencies
                r = AP.findResponses(lat,100);
                rates = [rates;r(:,1)];
            end
            
            null_mean = mean(rates);
            null_std = std(rates);
        end
        
        function [responses] = filter(AP, data, key, val)
            % filter responses by orientation, coefficient, type
            if strcmp(key,'orientation')
                type = 2;
            elseif strcmp(key,'type')
                type = 3;
            elseif strcmp(key,'coefficient')
                type = 4;
            elseif strcmp(key,'stroke')
                type = 6;
            end
            ind = find(data(:,type)==val);
            responses = data(ind,:);
        end
        
        function [responses] = getParabolaResponses(AP)
            responses = AP.filter(AP.responses,'type',0);
        end
        
        function [responses] = getAbsResponses(AP)
            responses = AP.filter(AP.responses,'type',1);
        end
        
        function [out] = findVariance(AP,varargin)
            
            p = inputParser;
            addRequired(p,'AP');
            addParamValue(p,'collapse_stroke',1);
            addParamValue(p,'responses',1);
            parse(p,AP,varargin{:});
            
            if p.Results.responses~=1
                raw_responses = p.Results.responses;
            else
                raw_responses = AP.responses;
            end
            
            
            % separate stimuli info from firing rates
            % warning: this code may be broken since inum column was added to
            % responses matrix
            if p.Results.collapse_stroke==1 % average across strokes
                stimuli = raw_responses(:,2:5);
            else
                stimuli = raw_responses(:,2:6);
            end
            
            responses = raw_responses(:,1);
            
            % find unique and duplicate stimuli
            [stimuli_unique, uix] = unique(stimuli,'rows','first');
            responses_unique = responses(uix);
            dupix = setdiff(1:size(stimuli,1),uix); %indices of duplicate rows
            
            % setup output variables
            out = {};
            out.counts = zeros(size(stimuli_unique,1),1);
            out.stds = zeros(size(stimuli_unique,1),1);
            out.means = zeros(size(stimuli_unique,1),1);
            out.hi_stds = [];
            out.rates = cell(size(stimuli_unique,1),1);
            
            % loop over each unique stimuli, then scan set of duplicates
            % for copies
            for i=1:size(stimuli_unique,1)
                rates = [responses_unique(i)];
                
                for j=1:length(dupix)
                    if stimuli(dupix(j),:)==stimuli_unique(i,:)
                        rates = [rates; responses(dupix(j))];
                    end
                end
                
                out.counts(i) = length(rates);
                if length(rates)==1
                    out.stds(i)=-10;
                else
                    out.stds(i) = std(rates);
                end
                
                if(length(rates) >= 4)
                    out.hi_stds = [out.hi_stds;std(rates)];
                end
                out.means(i) = mean(rates);
                
                out.rates{i} = rates;
            end
            
            n_empty = AP.N_STIMULI - size(stimuli_unique,1);
            c = out.counts;
            c((length(out.counts)+1):(length(out.counts)+1+n_empty)) = 0;
            [out.counts_counts, out.counts_bins] = hist(c,0:max(out.counts));
            
%             pie(out.counts_counts);
%             legend(cellfun(@num2str,num2cell(out.counts_bins),'UniformOutput',0));
        end

        function explainable_r = explainableVariance(AP)
        % EXPLAINABLEVARIANCE measures correlation b/w odd and even trials
        % It returns the pearson coefficient EXPLAINABLE_R. To find the 
        % explainable variance this should be squared.
        %
        % If no output args specified then plots a histogram showing the
        % different variances calculated on different resamplings.
        %
        % NOTE: This function differs from RasterUtil.explainableVariance()
        % in that it collapses across stroke before averaging different
        % responses to get a kernel estimate. This leads to higher
        % explainable variances that are much closer in line with the
        % explained variances, leading to reasonable prediction scores.
        % Using RasterUtil.explainableVariance often results in prediction 
        % scores higher than 1, and in many cases as high as 3

            N_DRAWS = 25;
            N_RESPONSES = size(AP.responses,1);
            rs = zeros(N_DRAWS,1);
            for i=1:N_DRAWS
            	odd_ix = randsample(N_RESPONSES,round(N_RESPONSES/2));
                even_ix = randsample(N_RESPONSES,round(N_RESPONSES/2));
                odd_kernel = AP.findKernel(AP.responses(odd_ix,:));
                even_kernel = AP.findKernel(AP.responses(even_ix,:));
                rs(i) = corr(odd_kernel(:),even_kernel(:),'rows','complete');
            end
            explainable_r = mean(rs);
            
            if nargout < 1
                figure();
                hist(rs);
                title(sprintf('Explainable Variance (mean r=%.2f)',explainable_r));
            end
        end
        
        function kernel = findKernel(AP, responses, varargin)
            p = inputParser;
            addRequired(p,'AP');
            addRequired(p,'responses');
            addParamValue(p,'significant',0);
            addParamValue(p,'method','bootstrap_t_test');
            parse(p,AP,responses,varargin{:});
            
            % If significant==0, then don't perform randomization analysis
            % to assess significance of the computed kernel
            if p.Results.significant==0
                
                kernel = NaN([length(AP.orientations),length(AP.coefficients),length(AP.types)]);
                
                % Find unique stimuli, and then use indices to grab
                % responses to each unique stimulus and then average them
                % together to create the kernel
                [unique_stims,~,ic] = unique(responses(:,2:4),'rows');
                rates = zeros(size(unique_stims,1),1);
                for i=1:size(unique_stims,1)
                    matches = responses(ic==i,1);
                    rates(i) = mean(matches);
                    theta = find(AP.orientations==unique_stims(i,1));
                    coeff = find(AP.coefficients==unique_stims(i,3));
                    type = find(AP.types==unique_stims(i,2));
                    kernel(theta,coeff,type)=rates(i);
                end
                
            else
                % Default: Bootstrap t-test
                if strcmp(p.Results.method,'bootstrap_t_test')

                    % Use Jackknifed Mean Kernel as Starting Point
                    [estimated_kernel, ~, ~] = AP.findJackknifedKernel(responses);
                    [lower_ci, upper_ci, mean_rate] = AP.findBootstrappedNull(responses);

                    kernel = estimated_kernel;
                    kernel((kernel < upper_ci) & (kernel > lower_ci)) = mean_rate;

                % Find Confidence Intervals via Bootstrap or Jackknife
                else
                    if strcmp(p.Results.method,'jackknife')
                        [estimated_kernel, ~, estimated_se] = AP.findJackknifedKernel(responses);
                    elseif strcmp(p.Results.method,'bootstrap')
                        [estimated_kernel, ~, estimated_se] = AP.findBootstrappedKernel(responses);
                    end

                    kernel = estimated_kernel;
                    kernel((estimated_kernel - estimated_se) <= 0) = 0;
                end
            end
        end

        function [kernel_mean, kernel_std, kernel_se] = findJackknifedKernel(AP, responses)

            N_SAMPLES = 20;
            N_TOTAL = size(responses,1);
            N_DISCARDS = floor(N_TOTAL/N_SAMPLES);
            
            full_kernel = AP.findKernel(responses);
            pseudo_kernels = zeros(N_SAMPLES, size(full_kernel,1), size(full_kernel,2), size(full_kernel,3));
            
            for i=1:N_SAMPLES
                discarded_i = N_DISCARDS*(i-1)+1:N_DISCARDS*i;
                fractional_responses = responses(setdiff(1:N_TOTAL, discarded_i),:);
                fractional_kernel = AP.findKernel(fractional_responses);
                pseudo_kernels(i,:,:,:) = (N_SAMPLES * full_kernel) - ((N_SAMPLES -1) * fractional_kernel);
            end
            
            kernel_mean = squeeze(nanmean(pseudo_kernels,1));
            kernel_std = squeeze(nanstd(pseudo_kernels,1));
            kernel_se = kernel_std / sqrt(N_SAMPLES);
        end
        
        function [kernel_mean, kernel_std, kernel_se] = findBootstrappedKernel(AP, responses)

            N_DRAWS = 10;
            N_TOTAL = size(responses,1);
            
            full_kernel = AP.findKernel(responses);
            bootstrapped_kernels = zeros(N_DRAWS, size(full_kernel,1), size(full_kernel,2), size(full_kernel,3));
            
            for i=1:N_DRAWS
                bootstrapped_responses = responses(randsample(N_TOTAL,N_TOTAL,true),:);
                bootstrapped_kernels(i,:,:,:) = AP.findKernel(bootstrapped_responses);
            end
            
            kernel_mean = squeeze(nanmean(bootstrapped_kernels,1));
            kernel_std = squeeze(nanstd(bootstrapped_kernels,1));
            kernel_se = kernel_std;
       end
        
        function [lower_bound, upper_bound, null_mean, null_dist] = findBootstrappedNull(AP, responses)

            N_DRAWS = 50;
            N_TOTAL = size(responses,1);
            
            full_kernel = AP.findKernel(responses);
            bootstrapped_kernels = zeros(N_DRAWS, size(full_kernel,1), size(full_kernel,2), size(full_kernel,3));
            
            for i=1:N_DRAWS
                % Resample responses (shuffling entire rows, with replacement)
                bootstrapped_responses = responses(randsample(N_TOTAL,N_TOTAL,true),:);
                
                % Scramble firing rates (shuffling only the firing rates, with replacement)
                bootstrapped_responses(:,1) = responses(randsample(N_TOTAL,N_TOTAL,true),1);
                
                bootstrapped_kernels(i,:,:,:) = AP.findKernel(bootstrapped_responses);
            end
            
            % Unravel kernels to increase cardinality of distribution
            null_dist = bootstrapped_kernels(:);
            
            % Find null mean firing rate
            null_mean = nanmean(null_dist);
            
            % Calculate cumulative distribution function
            [cdf_p,cdf_rate] = ecdf(null_dist);
            
            % Upper and lower confidence interval bounds
            lower_bound = cdf_rate(find(cdf_p <= 0.05,1,'last'));
            upper_bound = cdf_rate(find(cdf_p >= 0.95,1,'first'));
        end
        
        function [r_final, actual_rates, fitted_rates, final_rates] = findModelKernel(AP, model)
                        
            image_info = AP.pf.rec(1).params.IMAGE_INFO;
            
            thin = min(AP.strokes);
            thick = max(AP.strokes);
            
            thin_kernel = zeros(length(AP.orientations),length(AP.coefficients),length(AP.types));
            thick_kernel = zeros(length(AP.orientations),length(AP.coefficients),length(AP.types));
            % final_kernel = zeros(length(AP.orientations),length(AP.coefficients),length(AP.types));
            
            for i=1:AP.N_STIMULI
                
                info = image_info{i};                 % switch from 0- to 1-based
                tokens = strsplit(info,[char(9) ' ']); % split by tab+space
                
                stimulus = model.stimulusLoader.getByImageNumber(AP.pf, i-1);
                
                response = model.stimulate(stimulus);
                
                % Parabola == 0, Abs == 1
                if strcmp(tokens(2),'parabola')
                    type = 0;
                else
                    type = 1;
                end

                coefficient = str2num(tokens{3});
                stroke = str2num(tokens{4});
                rotation = str2num(tokens{5});
                ori = 90 - rotation - 360;
                                
                % save kernel
                if stroke==thin
                   thin_kernel = AP.setKernel(thin_kernel, ori, coefficient, type, response.subunit_sum);
                else
                   thick_kernel = AP.setKernel(thick_kernel, ori, coefficient, type, response.subunit_sum);
                end
            end
            
            final_kernel = (thin_kernel + thick_kernel)/2;
            [~, ~, thin_rates, ~] = AP.kernel2vectors(thin_kernel);
            [~, ~, thick_rates, ~] = AP.kernel2vectors(thick_kernel);
            [~, ~, final_rates, ~] = AP.kernel2vectors(final_kernel);
            [~, ~, actual_rates, ~] = AP.kernel2vectors(AP.findKernel(AP.responses), 'include_nan',1);
            thin_rates = thin_rates(~isnan(actual_rates));
            thick_rates = thick_rates(~isnan(actual_rates));
            final_rates = final_rates(~isnan(actual_rates));
            actual_rates = actual_rates(~isnan(actual_rates));
            
            r_thin = corr(actual_rates, thin_rates);
            r_thick = corr(actual_rates, thick_rates);
            r_final = corr(actual_rates, final_rates);
            r_stroke = corr(thin_rates, thick_rates);
                        
            f = @(x,xdata)x(1)*xdata.^x(2)+x(3);
            options = optimoptions('lsqcurvefit','Display','iter','TolFun',1e-10);
            x = lsqcurvefit(f,[1 1 0],final_rates,actual_rates,[],[],options);
            fprintf('f(x) = %.3fx^%.3f %+0.3f\n',x(1),x(2),x(3));
            fitted_rates = f(x,final_rates);
            r_fitted = corr(actual_rates,fitted_rates);
            
%             figure();
%             
%             subplot(2,2,1);
%             scatter(thin_rates, actual_rates);
%             xlabel('thin');
%             ylabel('actual');
%             title(sprintf('r^2 = %0.3f',r_thin^2));
%             axis square;
%             
%             subplot(2,2,2);
%             scatter(thick_rates, actual_rates);
%             axis square;
%             xlabel('thick');
%             ylabel('actual');
%             title(sprintf('r^2 = %0.3f',r_thick^2));
%             
%             subplot(2,2,3);
%             scatter(final_rates, actual_rates);
%             axis square;
%             xlabel('mean of thick and thin');
%             ylabel('actual');
%             title(sprintf('r^2 = %0.3f',r_final^2));
%             
%             subplot(2,2,4);
% %             scatter(thin_rates, thick_rates);
% %             axis square;
% %             xlabel('thin');
% %             ylabel('thick');
% %             title(sprintf('r^2 = %0.3f',r_stroke^2));
%             scatter(fitted_rates, actual_rates);
%             axis square;
%             xlabel('fitted');
%             ylabel('actual');
%             title(sprintf('r^2 = %0.3f',r_fitted^2));
        end
        
        function [r_unfitted, rates_observed, rates_predicted, vm_predicted] = findModelKernel2(AP, model)
        % FINDMODELKERNEL2 compares a model and the responses (as opposed
        % to findModelKernel1 which compares a model against the kernel)
            
            % Stimulate Model with Each Stimuli
            APL = AnglePlayLoader(AP.pf);
            N_STIMULI = length(AP.pf.rec(1).params.IMAGE_INFO);
            stimuli_responses = zeros(N_STIMULI,1);
            for i=1:N_STIMULI
                stimulus = APL.getByImageNumber(AP.pf, i-1);
                response = model.stimulate(stimulus);
                stimuli_responses(i) = response.subunit_sum;
            end
            
            % Generate Observed and Predicted Response Vectors
            rates_observed = AP.responses(:,1);
            vm_predicted = zeros(size(rates_observed));
            for i=1:AP.responses           
                inum = AP.responses(i,8);
                vm_predicted(i) = stimuli_responses(inum + 1); % switch from 0- to 1-based indexing
            end
            
            % Fit a Nonlinearity
            f = @(x,xdata)x(1)*xdata.^x(2)+x(3);
            options = optimoptions('lsqcurvefit','Display','iter','TolFun',1e-10);
            x = lsqcurvefit(f,[1 2 0],vm_predicted,rates_observed,[],[],options);
            fprintf('f(x) = %.3fx^%.3f %+0.3f\n',x(1),x(2),x(3));
            rates_predicted = f(x,vm_predicted);
            
            r_unfitted = corr(rates_observed, vm_predicted);
            r_fitted = corr(rates_observed,rates_predicted);
            
            figure();
            
            subplot(1,2,1);
            scatter(vm_predicted, rates_observed);
            xlabel('V_m Predicted');
            ylabel('Rates Observed');
            title(sprintf('Linear Model (r^2 = %0.3f)',r_unfitted^2));
            axis square;
           
            subplot(1,2,2);
%             scatter(thin_rates, thick_rates);
%             axis square;
%             xlabel('thin');
%             ylabel('thick');
%             title(sprintf('r^2 = %0.3f',r_stroke^2));
            scatter(rates_predicted, rates_observed);
            axis square;
            xlabel('Rates Predicted');
            ylabel('Rates Observed');
            title(sprintf('Linear + AnglePlay NL Model (r^2 = %0.3f)',r_fitted^2));
        end

        function r_full = findKernelCorrelation(AP, crossValidationRatio, varargin)
        % FINDKERNELCORRELATION computes prediction score using kernel
        % Cross-validation is performed by using a fraction of the data
        % (provided as the CROSSVALIDATIONRATIO, e.g. 0.9) as a training set
        % and then using the kernel computed from this training set to
        % predict responses to the remaining responses.
        %
        % If CROSSVALIDATIONRATIO==1, the function measures the
        % prediction of the full kernel against itself.
        %
        % If no output arguments are provided, scatter plots are generated
        %
        % [r_full, rates_observed] = findKernelCorrelation(AP, crossValidationRatio, varargin)
        %
        %  INPUT
        %    crossValidationRatio - % of data to use as training set [0,1]
        %
        %    options - specified as <'option', value> pairs:
        %      bootstrap - get better estimate of r_full via boostrapping
        %      sig - (BROKEN) test prediction score of significant kernel only
        %
        %  OUTPUT
        %    r_full - correlation between kernel prediction and observed
        %    rates
        %    (figure) if no output variables assigned
        
            p = inputParser;
            addRequired(p,'AP');
            addRequired(p,'crossValidationRatio');
            addParameter(p,'bootstrap',0);
            addParameter(p,'sig',0); % test prediction score of significant kernel only
            parse(p,AP,crossValidationRatio,varargin{:});
            
            % If bootstrap option provided then do analysis with many
            % resamplings to better estimate r
            if p.Results.bootstrap==1
                N_DRAWS = 100;
                scores = zeros(N_DRAWS,1);
                for i=1:N_DRAWS
                    scores(i) = AP.findKernelCorrelation(crossValidationRatio);
                end
                r_full = mean(scores);
                
                % Plot histogram of r
                if nargout < 1
                    figure();
                    hist(scores);
                    title(sprintf('Kernel Performance, Bootstrapped Results (crossValidationRatio = %.2f, mean r=%.2f)',crossValidationRatio,r_full));
                    xlabel('r');
                end

                return;
            end
        
            if crossValidationRatio==1
                training_data_i = 1:length(AP.responses);
                prediction_data_i = 1:length(AP.responses);
            else 
                n_responses = length(AP.responses);
                n_training_data = round(n_responses*crossValidationRatio);
                n_prediction_data = n_responses - n_training_data;
                training_data_i = randperm(length(AP.responses),n_training_data);
                prediction_data_i = setxor(1:length(AP.responses),training_data_i);
            end
            
            training_data = AP.responses(training_data_i,:);
            prediction_data = AP.responses(prediction_data_i,:);
            
            kernel_full = AP.findKernel(training_data,'significant',p.Results.sig);
%             kernel_sig = AP.findKernel(training_data,'significant',1);
            
            N_STIMULI = length(AP.pf.rec(1).params.IMAGE_INFO);
            
            % Generate Observed and Predicted Response Vectors
            rates_observed = prediction_data(:,1);
            rates_predicted_full = zeros(size(rates_observed));
%             rates_predicted_sig = zeros(size(rates_observed));
            for i=1:length(prediction_data)
                ori = prediction_data(i,2);
                coeff = prediction_data(i,4);
                type = prediction_data(i,3);
                rates_predicted_full(i) = AP.searchKernel(kernel_full,ori,coeff,type);
%                 rates_predicted_sig(i) = AP.searchKernel(kernel_sig,ori,coeff,type);
            end
            
            % Compute full prediction score!
            r_full = corr(rates_observed, rates_predicted_full,'rows','complete');
            
            % r_sig tests how well the 'significant' kernel predicts
            % responses (i.e. the kernel stripped of all entries that
            % that could be explained by null hypothesis)
%             r_sig = corr(rates_observed, rates_predicted_sig,'rows','complete');
%             autocorrelation = corr(rates_predicted_full, rates_predicted_sig,'rows','complete');
            
            if nargout < 1
                figure();

%                 subplot(1,3,1);
                scatter(rates_predicted_full, rates_observed);
                axis square;
                xlabel('Rates Predicted');
                ylabel('Rates Observed');
                title(sprintf('AnglePlay Kernel Performance (r^2 = %0.3f)',r_full^2));

%                 subplot(1,3,2);
%                 scatter(rates_predicted_sig, rates_observed);
%                 axis square;
%                 xlabel('Rates Predicted');
%                 ylabel('Rates Observed');
%                 title(sprintf('Significant AnglePlay Kernel (r^2 = %0.3f)',r_sig^2));
% 
%                 subplot(1,3,3);
%                 scatter(rates_predicted_sig, rates_predicted_full);
%                 axis square;
%                 xlabel('Significant');
%                 ylabel('Full');
%                 title(sprintf('Full vs Significant Kernel (r^2 = %0.3f)',autocorrelation^2));
            end
        end
        
        function [] = plotNullPDF(AP)
            
            [~,~,~,null_dist] = AP.findBootstrappedNull(AP.responses);
            
            kernel = AP.findKernel(AP.responses,'significant',0);
            kernel_par = kernel(:,:,1);
            dist_par = kernel_par(:);
            kernel_abs = kernel(:,:,2);
            dist_abs = kernel_abs(:);
            
            [y_null,x_null] = hist(null_dist,10);
            dx_null = diff(x_null(1:2));
            pdf_null = y_null / sum(y_null*dx_null);
            
            [y_par,x_par] = hist(dist_par,10);
            dx_par = diff(x_par(1:2));
            pdf_par = y_par / sum(y_par*dx_par);
            
            [y_abs,x_abs] = hist(dist_abs,10);
            dx_abs = diff(x_abs(1:2));
            pdf_abs = y_abs / sum(y_abs*dx_abs);
            
                        
%            figure();
%             plot(x_null,pdf_null,'--','Color',[0.7 0.7 0.7]);
%             hold on;
%             plot(x_par,pdf_par,'-r');
%             plot(x_abs,pdf_abs,'-b');
%             title('PDF Comparison');
%             legend('null','parabola','abs');
            
           % dist_abs = dist_par;
            
            
            [abs_cdf,abs_x2] = ecdf(dist_abs);
            [par_cdf,par_x2] = ecdf(dist_par);
            [null_cdf,null_x2] = ecdf(null_dist);
            
            figure();
            subplot(2,2,1);
            plot(abs_x2,abs_cdf,'-r');
            hold on;
            plot(par_x2,par_cdf,'-b');
            plot(null_x2,null_cdf,'--k');
            legend('abs','par','null');
            
            abs_cdf_s = smooth(abs_cdf,11,'rloess');
            
            subplot(2,2,2);
            plot(abs_x2, abs_cdf_s,'-o');
            x_max = max(abs_x2);
            title('smoothed');
            
            if(abs_x2(1)==0 && abs_x2(2)==0)
                abs_x2 = abs_x2(2:length(abs_x2));
                abs_cdf_s = abs_cdf_s(2:length(abs_cdf_s));
            end
            
            x3 = 0:5:x_max;

            cdf_int = interp1(abs_x2,abs_cdf_s,x3);
            
            subplot(2,2,3);
%             plot(0:1:x_max,cdf_int);
%             title('interpolated cdf');
            plot(x3,cdf_int);
            
            subplot(2,2,4);
            
            plot(x3(1:length(x3)-1), (diff(cdf_int)./diff(x3)));
            
            SMOOTHING = 15;
            DX = 10;
            
            [abs_cdf,abs_x2] = ecdf(dist_abs);
            abs_cdf_s = smooth(abs_cdf,SMOOTHING,'rloess');
            x_max = max(abs_x2);
            if(abs_x2(1)==0 && abs_x2(2)==0)
                abs_x2 = abs_x2(2:length(abs_x2));
                abs_cdf_s = abs_cdf_s(2:length(abs_cdf_s));
            end
            x3 = 0:DX:x_max;
            cdf_int = interp1(abs_x2,abs_cdf_s,x3);
            
            pdf_x = x3(1:length(x3)-1);
            pdf_y = (diff(cdf_int)./diff(x3));
            
            plot(pdf_x, pdf_y,'-r');
            hold on;
            
            
            dist_abs = null_dist;
            [abs_cdf,abs_x2] = ecdf(dist_abs);
            abs_cdf_s = smooth(abs_cdf,SMOOTHING,'rloess');
            x_max = max(abs_x2);
            if(abs_x2(1)==0 && abs_x2(2)==0)
                abs_x2 = abs_x2(2:length(abs_x2));
                abs_cdf_s = abs_cdf_s(2:length(abs_cdf_s));
            end
            x3 = 0:DX:x_max;
            cdf_int = interp1(abs_x2,abs_cdf_s,x3);
            
            pdf_x = x3(1:length(x3)-1);
            pdf_y = (diff(cdf_int)./diff(x3));
            
            plot(pdf_x, pdf_y,'--','Color',[0.7 0.7 0.7]);
            hold on;
            
            dist_abs = dist_par;
            [abs_cdf,abs_x2] = ecdf(dist_abs);
            abs_cdf_s = smooth(abs_cdf,SMOOTHING,'rloess');
            x_max = max(abs_x2);
            if(abs_x2(1)==0 && abs_x2(2)==0)
                abs_x2 = abs_x2(2:length(abs_x2));
                abs_cdf_s = abs_cdf_s(2:length(abs_cdf_s));
            end
            x3 = 0:DX:x_max;
            cdf_int = interp1(abs_x2,abs_cdf_s,x3);
            
            pdf_x = x3(1:length(x3)-1);
            pdf_y = (diff(cdf_int)./diff(x3));
            
            plot(pdf_x, pdf_y,'-b');
            hold on;
            
            legend('abs','null','par');
            
        end
        
        function rate = searchKernel(AP, kernel, ori, coeff, type)
            % Find indices
            ori_i = find(AP.orientations==ori);
            coeff_i = find(AP.coefficients==coeff);
            type_i = find(AP.types==type);
            
            % Lookup in the matrix
            rate = kernel(ori_i, coeff_i, type_i);
        end
        
        function kernel = setKernel(AP, kernel, ori, coeff, type, rate)
            % Find indices
            ori_i = find(AP.orientations==ori);
            coeff_i = find(AP.coefficients==coeff);
            type_i = find(AP.types==type);
            
            % Lookup in the matrix
            kernel(ori_i, coeff_i, type_i) = rate;
        end
        
        function [oris,coeffs,rates,types] = kernel2vectors(AP, kernel, varargin)
            
            p = inputParser;
            addRequired(p,'AP');
            addRequired(p,'kernel');
            addParamValue(p,'include_nan',0);
            parse(p,AP,kernel,varargin{:});
            
            oris = [];
            coeffs = [];
            rates = [];
            types = [];
            for i=1:length(AP.orientations)
                for j=1:length(AP.coefficients)
                    for k=1:length(AP.types)
                        rate = kernel(i,j,k);
                        if ~isnan(rate) || p.Results.include_nan==1
                            oris = [oris; AP.orientations(i)];
                            coeffs = [coeffs; AP.coefficients(j)];
                            rates = [rates; rate];
                            types = [types; AP.types(k)];
                        end
                    end
                end
            end
        end
                
        function [] = plotJackknifedTuning(AP, type)
            
            if strcmp(type,'abs')
                type = 1;
                full_responses = AP.getAbsResponses();
            else
                type = 0;
                full_responses = AP.getParabolaResponses();
            end
         
            full_kernel = AP.findKernel(full_responses);
            [fractional_kernels_mean, ~, fractional_kernels_se] = AP.findJackknifedKernel(full_responses);
            
            subplot(2,2,1);
            AP.plotTuning(full_kernel,type,'angles',0);
            title('Full Kernel','fontweight','bold');
                        
            subplot(2,2,2);
            AP.plotTuning(fractional_kernels_mean,type,'angles',0);
            title('Jackknifed Kernel Mean','fontweight','bold');
            colorbar();

            subplot(2,2,3);
            AP.plotTuning(fractional_kernels_se,type,'angles',0);
            title('Jackknifed Kernel SE','fontweight','bold');
            colorbar();
            
            subplot(2,2,4);
            
            sig_kernel = full_kernel;
            sig_kernel((full_kernel - fractional_kernels_se) <= 0) = 0;
            AP.plotTuning(sig_kernel,type,'angles',0);
            title('Significant Kernel','fontweight','bold');
            
%             difference = fractional_kernels_mean - full_kernel;
%             AP.plotTuning(difference,type,'angles',0);
%             title('Jackknifed - Full Kernel','fontweight','bold');
%             colorbar();
%             diff_range = max(abs(min(difference(:))),abs(max(difference(:))));
%             caxis([-1*diff_range diff_range]);
            
            boxtitle([AP.plot_title ' Jackknife Analysis']);
        end
        
        function [] = plotRandomizedTuning(AP,type)
            figure();
             if strcmp(type,'abs')
                type = 1;
                full_responses = AP.getAbsResponses();
            else
                type = 0;
                full_responses = AP.getParabolaResponses();
            end
         
            full_kernel = AP.findKernel(full_responses);
            [fractional_kernels_mean, ~, fractional_kernels_se] = AP.findJackknifedKernel(full_responses);
            
            subplot(2,2,1);
            AP.plotTuning(full_kernel,type,'angles',0);
            title('Full Kernel','fontweight','bold');
                        
            subplot(2,2,2);
            sig_kernel = fractional_kernels_mean;
            sig_kernel((fractional_kernels_mean - fractional_kernels_se) <= 0) = 0;
            AP.plotTuning(sig_kernel,type,'angles',0);
            title('Post-Jackknife Kernel','fontweight','bold');
            
            subplot(2,2,3);
            [lower_ci, upper_ci, null_mean] = AP.findBootstrappedNull(full_responses);
            
            final_kernel = sig_kernel;
            final_kernel((final_kernel < upper_ci) & (final_kernel > lower_ci)) = 0;
            AP.plotTuning(final_kernel,type,'angles',0);
            title('T-Tested Kernel','fontweight','bold');
            
            boxtitle([AP.plot_title ' Randomization Analysis']);
            
        end
        
        function [] = plotBootstrappedTuning(AP, type)
            
            if strcmp(type,'abs')
                type = 1;
                full_responses = AP.getAbsResponses();
            else
                type = 0;
                full_responses = AP.getParabolaResponses();
            end
         
            full_kernel = AP.findKernel(full_responses);
            [fractional_kernels_mean, ~, fractional_kernels_se] = AP.findBootstrappedKernel(full_responses);
            
            subplot(2,2,1);
            AP.plotTuning(full_kernel,type,'angles',0);
            title('Full Kernel','fontweight','bold');
                        
            subplot(2,2,2);
            AP.plotTuning(fractional_kernels_mean,type,'angles',0);
            title('Bootstrapped Kernel Mean','fontweight','bold');
            colorbar();
            
            subplot(2,2,3);
            AP.plotTuning(fractional_kernels_se,type,'angles',0);
            title('Bootstrapped Kernel SE','fontweight','bold');
            colorbar();
            
            subplot(2,2,4);
            difference = fractional_kernels_mean - full_kernel;
            AP.plotTuning(difference,type,'angles',0);
            title('Bootstrapped - Full Kernel','fontweight','bold');
            colorbar();
            diff_range = max(abs(min(difference(:))),abs(max(difference(:))));
            caxis([-1*diff_range diff_range]);
            
            boxtitle([AP.plot_title ' Bootstrap Analysis']);
        end
        
        function [xs,ys,zs] = plotParabolaTuning(AP,varargin)
            
            p = inputParser;
            addRequired(p,'AP');
            addParamValue(p,'significant',1);
            addParamValue(p,'unit',1);
            addParamValue(p,'zscore',1);
            parse(p,AP,varargin{:});
            
            if strcmp(p.Results.unit,'angles')
                unit = 'angles';
            else
                unit = 'curvature';
            end
            
          %  responses = AP.getParabolaResponses();
            kernel = AP.findKernel(AP.responses,'significant',p.Results.significant);
            
            [xs,ys,zs] = AP.plotTuning(kernel,0,unit,p.Results.zscore);
        end
        
        function [xs,ys,zs] = plotAbsTuning(AP,varargin)
            
            p = inputParser;
            addRequired(p,'AP');
            addParamValue(p,'significant',1);
            addParamValue(p,'unit',1);
            addParamValue(p,'zscore',1); % (0/1)
            parse(p,AP,varargin{:});
            
            if strcmp(p.Results.unit,'curvature')
                unit = 'curvature';
            else
                unit = 'angles';
            end
                
            % instead of plotting coefficients, plot using angles
           % responses = AP.getAbsResponses();
            kernel = AP.findKernel(AP.responses,'significant',p.Results.significant);
            
            [xs,ys,zs] = AP.plotTuning(kernel,1,unit,p.Results.zscore);
        end
        
        function [xs,ys,zs] = plotTuning(AP,kernel,data_type,unit,zscore)
            
            % ignore kernel entries for other type
            if data_type==0 %parabola
                kernel(:,:,2)=NaN;
            else %abs
                kernel(:,:,1)=NaN;
            end
            [oris, coeffs, rates, types] = AP.kernel2vectors(kernel);
            
            % Z-Score the Firing Rates Using Null Data
            if isempty(zscore) || zscore==1
                max_rate = (AP.max_rate - AP.null_mean)/AP.null_std;
                rates = (rates - AP.null_mean)/AP.null_std;
                axis_range = [1 max_rate];
                
            % Use Firing Rates (s/s)
            else
                max_rate = max(rates(:));
                axis_range = [0 max_rate];
                %mean_rate = mode(rates(:));
                %axis_range = [(mean_rate - (max_rate - mean_rate)) max_rate];
            end
            
            if strcmp(unit,'angles')
                
                if data_type==1 % abs
                     angles = 2*atand(1./AP.coefficients); % angles of entire set 
                     angs = 2*atand(1./coeffs); % angles from response set
                else % parabola: find where they intercept a circle with r=2
                    syms x positive;
                    xs = zeros(length(coeffs),1);
                    for i=1:length(coeffs)
                        xs(i) = double(solve(coeffs(i)*x^2 == sqrt(-x^2 + 4),x));
                    end
                    ys = coeffs.*(xs.^2);
                    angs = 2*atand(xs./ys);
                    clear x;
                end
                              
                
                [xs, ys] = meshgrid(min(AP.orientations):1:max(AP.orientations),0:1:180);
                zs = griddata(oris,angs,rates,xs,ys,'cubic');

                h = contourf(xs,ys,zs,50);
                try
                    shading flat;
                catch err
                end
                hold on;
                scatter(oris,angs,'.','MarkerEdgeColor','w');%[0.5 0.5 0.5]);
                
                ylabel('Angle');
                set(gca,'YDir','reverse');
                
            else
                curvatures = AP.coefficients * 2;
                curvs = coeffs * 2;

                [xs, ys] = meshgrid(min(AP.orientations):0.5:max(AP.orientations),min(curvatures):0.01:max(curvatures));

                zs = griddata(oris,curvs,rates,xs,ys,'cubic');

                contourf(xs,ys,zs,50);
                try
                    shading flat;
                catch err
                end
                hold on;
                scatter(oris,curvs,'+','MarkerEdgeColor','w');%[0.5 0.5 0.5]);
                
                ylabel('Curvature');
                
            end
            %colormap(hotcold(100));
            caxis(axis_range);
            xlabel('Orientation');
            colorbar();
        end
        
        function [sig_par_rates,sig_abs_rates] = plotTypes(AP)
            responses = AP.responses;
            
            oris = [];
            coeffs = [];
            rates = [];
            par_rates = [];
            abs_rates = [];
            for theta=1:length(AP.orientations)
                ori_matches = AP.filter(responses,'orientation',AP.orientations(theta));
                for a=1:length(AP.coefficients)
                    coeff_matches = AP.filter(ori_matches,'coefficient',AP.coefficients(a));
                    
                    abs_matches = AP.filter(coeff_matches,'type',1);
                    par_matches = AP.filter(coeff_matches,'type',0);
                    
                    par_rate = mean(par_matches(:,1));
                    abs_rate = mean(abs_matches(:,1));
                    if ~isnan(par_rate) && ~isnan(abs_rate)
                        par_rates = [par_rates; par_rate];
                        abs_rates = [abs_rates; abs_rate];
                    end
                end
            end
            
            %% Plot Each Dot
            % if significant, plot as blue else as gray
            sig_par_rates = []; % significant parabola rates
            sig_abs_rates = []; % significant abs rates
            for i=1:length(par_rates)
                if((par_rates(i) > AP.null_mean) || (abs_rates(i) > AP.null_mean))
                    marker_color = 'blue';
                    sig_par_rates = [sig_par_rates; par_rates(i)];
                    sig_abs_rates = [sig_abs_rates; abs_rates(i)];
                else
                    marker_color = [0.5 0.5 0.5];
                end
                scatter(par_rates(i),abs_rates(i),'MarkerEdgeColor',marker_color);
                hold on;
            end
            
            %% Stats
            rate_max = max(max(par_rates(:)),max(abs_rates(:)));
            
            % correlation of all points
            r_all = corr(par_rates,abs_rates);
            
            % correlation of only points with rates higher than acausal mean
            r_significant = corr(sig_par_rates,sig_abs_rates); 
            
            
            %% Finish Plot
            
            % Plot Acausal Mean Boundary
            line([AP.null_mean AP.null_mean],[0 AP.null_mean],'LineStyle',':','Color',[ 0.5 0.5 0.5]);
            line([0 AP.null_mean],[AP.null_mean AP.null_mean],'LineStyle',':','Color',[ 0.5 0.5 0.5]);

            % Plot Diagonal Line ~ r=1
            line([0 rate_max],[0 rate_max],'Color','k');
  
            xlabel('Parabola (s/s)');
            ylabel('Absolute Value (s/s)');
            title(sprintf('Parabola vs Abs (r=%.2f, r_{sig}=%.2f)',r_all,r_significant),'fontweight','bold');
            axis([0 rate_max 0 rate_max]); axis square;
        end
        
        function [] = plotStrokes(AP)
            responses = AP.responses;
            
            THIN = min(AP.strokes);
            THICK = max(AP.strokes);
            
            oris = [];
            coeffs = [];
            rates = [];
            thin_rates = [];
            thick_rates = [];
            for theta=1:length(AP.orientations)
                ori_matches = AP.filter(responses,'orientation',AP.orientations(theta));
                for a=1:length(AP.coefficients)
                    coeff_matches = AP.filter(ori_matches,'coefficient',AP.coefficients(a));
                    
                    thin_matches = AP.filter(coeff_matches,'stroke',THIN);
                    thick_matches = AP.filter(coeff_matches,'stroke',THICK);
                    
                    thin_rate = mean(thick_matches(:,1));
                    thick_rate = mean(thin_matches(:,1));
                    if ~isnan(thin_rate) && ~isnan(thick_rate)
                        thin_rates = [thin_rates; thin_rate];
                        thick_rates = [thick_rates; thick_rate];
                    end
                end
            end

            scatter(thin_rates,thick_rates);
            rate_max = max(max(thin_rates(:)),max(thick_rate(:)));
           
            hold on;
            line([0 rate_max],[0 rate_max],'Color','k');
            xlabel(['Thin (stroke=' num2str(THIN) ')']);
            ylabel(['Thick (stroke=' num2str(THICK) ')']);
            title(['Stroke Responses (r=' num2str(corr(thin_rates,thick_rates)) ')'],'fontweight','bold');
            axis([0 rate_max 0 rate_max]); axis square;
        end
        
        function [] = plotCrossValidation(AP)
            
            N_DRAWS = 20;
            N_TOTAL = size(AP.responses,1);
            
            % Generate Fractional P-Thetas
            fractions = 0.2:0.1:1;
            
            correlations = zeros(N_DRAWS,length(fractions),1);
            sig_correlations = zeros(N_DRAWS,length(fractions),1); % sig == significant
            
            % Measure Ground Truth (aka full-dataset)
            full_responses = AP.responses;
            true_kernel = AP.findKernel(full_responses);
            
            
            
            for i=1:length(fractions)
                N_FRACTIONAL_RESPONSES = floor(N_TOTAL*fractions(i));
                
                for draw=1:N_DRAWS
                    % Generate fractional kernel
                    rand_indices = randperm(N_TOTAL,N_FRACTIONAL_RESPONSES);
                    fractional_responses = full_responses(rand_indices,:);
                    fractional_kernel = AP.findKernel(fractional_responses);
                    
                    % Find indices of kernel entries that are not NaN in
                    % both the full and fractional vectors
                    match_i = find((~isnan(true_kernel(:)) & ~isnan(fractional_kernel(:)))==1);
                    true_matches = true_kernel(match_i);
                    fractional_matches = fractional_kernel(match_i);
                    sig_true_matches = true_matches(find(true_matches >= AP.null_mean));
                    sig_fractional_matches = fractional_matches(find(true_matches >= AP.null_mean));

                    % Compute Correlations
                    correlations(draw,i) = corr(true_matches, fractional_matches);
                    sig_correlations(draw,i) = corr(sig_true_matches,sig_fractional_matches);
                end
            end
       
            % convert fractions to n trials
            n_trials = round(fractions*size(AP.pf.rec,2));
            
            % convert r to r^2
            r_squared = mean(correlations,1).^2;
            sig_r_squared = mean(sig_correlations,1).^2;

            plot(n_trials,r_squared,'-o');
            hold on;
            plot(n_trials,sig_r_squared,'--ro');
            title('Cross-Validation (r^2)','fontweight','bold');
            xlim([0 max(n_trials)]);
            xlabel('trials');
            ylabel('r^2');
            legend('all rates','rates > acausal mean','Location','SouthEast');
        end
        
        function [] = loadBlurredStimulusImages(AP,kernel,imagedir)
            SLICE_FRACTION = 0.8;
            RAW_IMG_SIZE = 800;
            KERNEL_STROKE_RATIO = 2;
            [oris,coeffs,~,types] = AP.kernel2vectors(kernel,'include_nan',1);
            
            calibration_img = AP.loadStimulusImage(imagedir, coeffs(1), types(1));
            IMG_SIZE = max(size(calibration_img)); % size of image loaded from file
            AP.stimulusSize = IMG_SIZE;
            AP.slice_fraction = SLICE_FRACTION;
            SLICE_SIZE = round(IMG_SIZE * SLICE_FRACTION); % used to eliminate edge effects
            
            % Create blur kernel that is the same width as the stroke used
            % in the images. 
            
            % stroke width as fraction of image size. uses raw image size
            % since the stroke width defined in the params file is rendered
            % in the originals and then scaled in other files
            STROKE_IMG_RATIO = max(AP.strokes)/RAW_IMG_SIZE;
            KERNEL_SIZE = KERNEL_STROKE_RATIO * STROKE_IMG_RATIO * IMG_SIZE; % in pixels
            blur_kernel = fspecial('disk',KERNEL_SIZE);
            
            AP.blurredStimulusImages = cell(length(oris),1);

            for i=1:length(oris)
                im = AP.loadStimulusImage(imagedir, coeffs(i), types(i));
                im_inverted = imcomplement(im);
                im_rotated = imrotate(im_inverted, oris(i));
                im_blurred = imfilter(im_rotated,blur_kernel,'replicate');
                rotated_dim = size(im_rotated,1);
                min_pad = round((rotated_dim - SLICE_SIZE)/2);
                im_cropped = im_blurred(min_pad:(min_pad + SLICE_SIZE - 1),min_pad:(min_pad + SLICE_SIZE - 1));
                
                % correct for differences in image intensity by normalizing
                % sum of all pixels to 1
                im_n = double(im_cropped)/sum(im_cropped(:));
                AP.blurredStimulusImages{i} = im_n;
            end
        end
        
        function file_n = findImageFileNumber(AP, coeff, type)
        % FINDIMAGEFILENUMBER finds the number of the image file, and is
        % used to load the image path
        
            r = sortrows(AP.responses,7);
            % loads thickest image
            for i=1:size(r,1)
                if (r(i,4)==coeff && r(i,3)==type)
                    file_n = r(i,7);
                end
            end
        end

        function pic = loadStimulusImage(AP, directory, coeff, type)
            img_file_n = AP.findImageFileNumber(coeff, type);
            filename = sprintf('image-%06d.png',img_file_n);
            path = [directory '/' filename];
            pic = imread(path,'png');
        end

        function [stim_avg, sta, final, rates_zscored] = plotImageDomain(AP,varargin)
        % PLOTIMAGEDOMAIN visualizes results as the sum of blurred and weighted
        % stimulus images.
        %    OPTION significant     (int) 0/1 show only significant values
        %    OPTION show_plots      (int) 0/1 show plots?
        %    OPTION show_null       (int) 0/1 show stimulus average plot?
        %    OPTION show_sta        (int) 0/1 show sta plot?
        %    OPTION show_final      (int) 0/1 show sta / stim avg plot?
        %    OPTION plot_parabolas  (int) 0/1 include parabolas in plot?
        %    OPTION plot_abs        (int) 0/1 inclue abs in plot?
        %    OPTION imagedir        (string) directory of stimulus images
        %    RETURN acc             (matrix) accumulator / result image
            
            p = inputParser;
            addRequired(p,'AP');
            addParamValue(p,'significant',1);
            addParamValue(p,'imagedir',AP.imagedir);
            addParamValue(p,'show_plots',1);
            addParamValue(p,'show_null',0);
            addParamValue(p,'show_sta',0);
            addParamValue(p,'show_final',0);
            addParamValue(p,'plot_parabolas',1);% 0 or 1
            addParamValue(p,'plot_abs',1);% 0 or 1
            parse(p,AP,varargin{:});
            
            SHOW_PLOTS = p.Results.show_plots;
            SHOW_NULL = p.Results.show_null;
            SHOW_STA = p.Results.show_sta;
            SHOW_FINAL = p.Results.show_final;
            
            if(SHOW_NULL==0 && SHOW_STA==0 && SHOW_FINAL==0)
                SHOW_FINAL = 1;
            end
        
            PLOT_PARABOLAS = p.Results.plot_parabolas;
            PLOT_ABS = p.Results.plot_abs; 
            ZSCORE_THRESHOLD = -10;
            ZSCORE_RECTIFY = 1;
            COLOR_PALETTE = 'jet'; % ('jet' or 'gray')
            

            kernel = AP.findKernel(AP.responses,'significant',p.Results.significant);
            [oris,coeffs,rates,types] = AP.kernel2vectors(kernel,'include_nan',1);
            
            % Load Blurred Stimulus Images
            if isempty(AP.blurredStimulusImages)
                AP.loadBlurredStimulusImages(kernel,p.Results.imagedir);
            end
            
            SLICE_SIZE = max(size(AP.blurredStimulusImages{1}));
            
            % z-score rates
            rates_zscored = (rates - AP.null_mean)/AP.null_std;
            if ZSCORE_RECTIFY==1
                rates_zscored(rates_zscored < 0) = 0; % remove negative values
            end
            max_z = max(rates_zscored);
            
            
           
            %% Stimulus Average
            
            stim_avg_accumulator = zeros(SLICE_SIZE,'double');

            for i=1:length(oris)
                
                if PLOT_PARABOLAS==0 && types(i)==0
                    continue;
                elseif PLOT_ABS==0 && types(i)==1
                    continue;
                end
                stim_avg_accumulator = stim_avg_accumulator + AP.blurredStimulusImages{i};
            end
            
            % Normalize so that the sum == 1
            stim_avg = stim_avg_accumulator / sum(stim_avg_accumulator(:));

            %% STA
            sta_accumulator = zeros(SLICE_SIZE,'double');
                        
            for i=1:length(oris)
                
                if isnan(rates_zscored(i))
                    continue;
                elseif PLOT_PARABOLAS==0 && types(i)==0
                    rates_zscored(i)=0;
                    continue;
                elseif PLOT_ABS==0 && types(i)==1
                    rates_zscored(i)=0;
                    continue;
                end
                
                im_scaled = AP.blurredStimulusImages{i} .* rates_zscored(i);
                if abs(rates_zscored(i)) > ZSCORE_THRESHOLD
                    sta_accumulator = sta_accumulator + im_scaled;
               	end
            end
            
            % Normalize so that the sum == 1
            sta = sta_accumulator / sum(sta_accumulator(:));
            
            %% Calculate Contrast
            % (a - b) / (a + b)
            
            % final = sta ./ stim_avg;
            final = (sta - stim_avg) ./ (sta + stim_avg);
            final(isnan(final)) = 0;
            
            % params for 

            %% Plot Results
            if SHOW_PLOTS==1
                
                % calculate padding used to get rid of dead zones caused by
                % image rotation
                padding = (AP.stimulusSize - (AP.stimulusSize * AP.slice_fraction)) / 2;
                
                % calculate position of rounded rectangle showing RF sigma
                if AP.task_version < 938
                    rf_pos = [(AP.stimulusSize/6) - padding, (AP.stimulusSize/6)-padding, ((2*AP.stimulusSize)/3), ((2*AP.stimulusSize)/3)];
                    rf_color = [1,1,1];
                    rf_style = ':';
                    rf_stroke = 2;
                else
                    rf_sigma = (AP.stimulusSize/10);
                    rf_width = 2 * rf_sigma;
                    rf_pos = [(AP.stimulusSize/2) - (rf_width/2) - padding, (AP.stimulusSize/2) - (rf_width/2) - padding, rf_width, rf_width];
                    rf_color = [0,0,0];
                    rf_style = '-';
                    rf_stroke = 1;
                    
                    % Calculate Alpha Envelope
                    nsigma = AP.params.nsigma; % @todo implement nsigma
                end
                
                if SHOW_NULL == 1
                    imagesc(stim_avg);
                    hold on;
                    rectangle('Position',rf_pos,'Curvature',[1 1],'LineWidth',rf_stroke,'EdgeColor',rf_color,'LineStyle',rf_style);
                    title('Stimulus Avg','fontweight','bold');
                    axis square;
                    set(gca,'XTick',[], 'YTick',[]);
                    colorbar();
                end

                if SHOW_STA == 1
                    imagesc(sta);
                    hold on;
                    rectangle('Position',rf_pos,'Curvature',[1 1],'LineWidth',rf_stroke,'EdgeColor',rf_color,'LineStyle',rf_style);

                    title('STA','fontweight','bold');
                    axis square;
                    set(gca,'XTick',[], 'YTick',[]);
                    colorbar();
                end

                if SHOW_FINAL == 1
                    h = imagesc(final);
                    axis square;
                    
                    % Add Gaussian Envelope in the Alpha Channel
                    % Gaussian calculated w/ sigma = nsigma * rf_sigma
                    if AP.task_version >= 938
                        MIN_OPACITY = 0.01; % when == 0, plot fades to white
                        GAUSSIAN_MULTIPLIER = 1.5; % makes a broader region be fully visible
                        
                        cdata = getimage(h);
                        frame_size_withcolor = size(cdata);
                        frame_size = frame_size_withcolor(1:2); % strip out rgb axis
                        
                         % undo edge slicing to figure out how big the full
                         % stimulus would be if it were displayed in the
                         % little tiny frame (used to calculate sigma)
                        frameStimulusSize = frame_size / AP.slice_fraction;
                        frameSigma = max(frameStimulusSize) / 10;
                        frameEnvelopeSigma = round(nsigma * frameSigma);
                        gaussian = fspecial('Gaussian',frame_size,frameEnvelopeSigma);
                        gaussian_n = gaussian / max(gaussian(:));
                        gaussian_multiplied = gaussian_n * GAUSSIAN_MULTIPLIER;
                        gaussian_multiplied(gaussian_multiplied > 1) = 1;
                        
                        gaussian_scaled = round((1 - MIN_OPACITY) * ((-255*gaussian_multiplied)+256));
                        set(h,'AlphaData',gaussian_scaled);
                    end

                    hold on;
                    rectangle('Position',rf_pos,'Curvature',[1 1],'LineWidth',rf_stroke,'EdgeColor',[1 1 1],'LineStyle',rf_style);

                    title('STA (Contrast)','fontweight','bold');
                    
                    set(gca,'XTick',[], 'YTick',[]);
                    % colormap(hotcold(100));
                    caxis([-0.5 0.5]);
                    colorbar();
                end
            end
        end
        
        function [] = plotPixelDomain(AP)
            
            PLOT_PARABOLAS = 1; % 0 or 1
            PLOT_ABS = 1; % 0 or 1
            ZSCORE_THRESHOLD = 2;
            COLOR_PALETTE = 'jet'; % ('jet' or 'gray')
            
            x_grid = -3:0.05:3;
            y_grid = -3:0.05:3;
            
            abs_responses = AP.getAbsResponses();
            abs_kernel = AP.findKernel(abs_responses);
            [abs_oris, abs_coeffs, abs_rates] = AP.kernel2vectors(abs_kernel);
            abs_types = ones(length(abs_oris),1);

            par_responses = AP.getParabolaResponses();
            par_kernel = AP.findKernel(par_responses);
            [par_oris, par_coeffs, par_rates] = AP.kernel2vectors(par_kernel);
            par_types = zeros(length(par_oris),1);
            
            unsorted_oris = [abs_oris; par_oris];
            unsorted_coeffs = [abs_coeffs; par_coeffs];
            unsorted_types = [abs_types; par_types];
            unsorted_rates = [abs_rates; par_rates];
            
            % sort kernel by firing rate, ascending
            kernel = sortrows([unsorted_oris,unsorted_coeffs,unsorted_types,unsorted_rates],4);
            
            oris = kernel(:,1);
            coeffs = kernel(:,2);
            types = kernel(:,3);
            rates = kernel(:,4);
            
            N_STIMULI = size(kernel,1);
            
            % points is 3d array (n_frames * X-coordinates * Y-coordinates)
            points = cell(N_STIMULI,1);
           
            % Generate Points for the Plot
            for i=1:N_STIMULI
                if types(i)==0  % Parabola
                    ys = coeffs(i) * x_grid.^2;
                else % Abs
                    ys = coeffs(i) * abs(x_grid);
                end
                
                % Rotate points, save for later
                [new_xs, new_ys] = AP.rotatePoints(x_grid, ys, oris(i));
                points{i} = [new_xs; new_ys];
            end

            % z-score rates
            max_rate = (AP.max_rate - AP.null_mean)/AP.null_std;
            rates_zscored = (rates - AP.null_mean)/AP.null_std;
            rates_zscored(rates_zscored < 0) = 0; % remove negative values
            max_z = max(rates_zscored);
            
            % Set colormap
            if strcmp(COLOR_PALETTE,'gray')==1
                colors = sortrows(gray(100),-1);
            else
                colors = jet(100);
            end
            
            % Plot Each Stimuli
            for i=1:N_STIMULI
                
                if PLOT_PARABOLAS==0 && types(i)==0
                    continue;
                elseif PLOT_ABS==0 && types(i)==1
                    continue;
                end
                
                if rates_zscored(i) > ZSCORE_THRESHOLD
                    color = colors(max(1,round(100*(rates_zscored(i)/max_z))),:);
                    xy = points{i};
                    plot(xy(1,:), xy(2,:), 'Color', color, 'LineWidth', 2);
                    hold on;
                end
            end
            axis([-2 2 -2 2]);
            axis square;
            colormap(colors);
            caxis([0 max_z]);
            set(gca,'XTick',[], 'YTick',[]);
            colorbar();
        end
        
        function [] = plotHistory(AP)
            
            ZSCORE_THRESHOLD = 2;
            latencies = 0:10:200;
            winsize = 30;
           
            % get the sta first so we know how big the images we'll be
            % saving are
            [~,~,sta] = AP.plotImageDomain('show_plots',0);
            
            % data structures to hold stas
            full_stas = zeros(length(latencies),size(sta,1), size(sta,2));
            par_stas = zeros(size(full_stas));
            abs_stas = zeros(size(full_stas));
            
            % arrays to hold variances
            full_var = zeros(length(latencies),1);
            par_var = zeros(length(latencies),1);
            abs_var = zeros(length(latencies),1);
            
            % arrays to hold weights (e.g. rates_zscored sums)
            full_weights = zeros(length(latencies),1);
            par_weights = zeros(length(latencies),1);
            abs_weights = zeros(length(latencies),1);
            
            
            for i=1:length(latencies)
                AP.responses = AP.findResponses(latencies(i), winsize);
                
                % store stas
                [~,~,full_stas(i,:,:),full_rates] = AP.plotImageDomain('show_plots',0);
                [~,~,par_stas(i,:,:),par_rates] = AP.plotImageDomain('show_plots',0,'plot_abs',0);
                [~,~,abs_stas(i,:,:),abs_rates] = AP.plotImageDomain('show_plots',0,'plot_parabolas',0);
                
                % store variance
                full_var(i) = var(full_stas(i,:));
                par_var(i) = var(par_stas(i,:));
                abs_var(i) = var(abs_stas(i,:));
                
                % store weights
                full_weights(i) = sum(full_rates > ZSCORE_THRESHOLD);
                par_weights(i) = sum(par_rates > ZSCORE_THRESHOLD);
                abs_weights(i) = sum(abs_rates > ZSCORE_THRESHOLD);
                
            end
           
            % shift latencies by winsize/2, so that each point corresponds
            % to the midpoint of the bin
            shifted_latencies = latencies + round(winsize/2);
           
            % normalize weights
            psth = phist(AP.pf);
            psth_min_i = find(psth(:,1) == min(shifted_latencies));
            psth_max_i = find(psth(:,1) == max(shifted_latencies));
            psth_slice = psth(psth_min_i:psth_max_i,:);
            full_weights_interpolated = interp1(shifted_latencies,full_weights,psth_slice(:,1));
            par_weights_interpolated = interp1(shifted_latencies,par_weights,psth_slice(:,1));
            abs_weights_interpolated = interp1(shifted_latencies,abs_weights,psth_slice(:,1));
            
%             psth_slice_n = psth_slice(:,2) / max(psth_slice(:,2));
%             
%             full_weights_n = full_weights_interpolated ./ psth_slice_n;
%             max_full_weight = max(full_weights_n);
%             
%             % normalize again so that max is 1
%             full_weights_n = full_weights_n / max_full_weight;
%             par_weights_n = (par_weights_interpolated ./ psth_slice_n) / max_full_weight;
%             abs_weights_n = (abs_weights_interpolated ./ psth_slice_n) / max_full_weight;
                                   
            figure();
            subplot(3,1,1);
            phist(AP.pf);
            title('PSTH','fontweight','bold');
            xlim([-100 400]);
            
            subplot(3,1,2);
            plot(shifted_latencies,full_weights,'-k');
            hold on;
            plot(shifted_latencies,abs_weights,'-r');
            plot(shifted_latencies,par_weights,'-b');
            legend('All','Abs','Parabolas');
            title('Weights','fontweight','bold');
            ylabel('sum of zscored rates');
            xlim([-100 400]);
            
%             subplot(4,1,3);
%             plot(psth_slice(:,1),full_weights_n,'-k');
%             hold on;
%             plot(psth_slice(:,1),abs_weights_n,'-r');
%             plot(psth_slice(:,1),par_weights_n,'-b');
%             legend('All','Abs','Parabolas');
%             title('Normalized Weights','fontweight','bold');
%             ylabel('sum of zscored rates / psth');
%             axis([-100 400 0 1]);
            
            subplot(3,1,3);
            plot(shifted_latencies,full_var,'-k');
            hold on;
            plot(shifted_latencies,abs_var,'-r');
            plot(shifted_latencies,par_var,'-b');
            legend('All','Abs','Parabolas');
            title('Kernel Information Content','fontweight','bold');
            ylabel('sta variance');
            xlim([-100 400]);
            
            figure();
            
            for i=1:length(latencies)
                
                subplot(4,6,i);
                start = latencies(i);
                stop = latencies(i)+winsize;
                
%                 % multiply by normalized weight
%                 midpoint = (start+stop)/2;
%                 weight_i = find(psth_slice(:,1)==midpoint);
%                 sta_n = full_weights_n(weight_i) * squeeze(full_stas((2*i-1),:,:));
%                 imagesc(sta_n);
                
                imagesc(squeeze(full_stas(i,:,:)));
                caxis([-0.5 0.5]);
                title([num2str(start) 'ms - ' num2str(stop) 'ms'],'fontweight','bold');
                set(gca,'XTick',[],'YTick',[]);
                axis square;
            end
        end
        
    end
    
    methods(Static)
        
        function [new_xs, new_ys] = rotatePoints(xs,ys,theta)
        % ROTATEPOINTS rotates a set of points (XS and YS) *counterclockwise* 
        % around the origin by some degrees THETA
            
           V = [xs;ys];
           R = [[cosd(theta),-sind(theta)];[sind(theta),cosd(theta)]];
           V_r = R*V;
           new_xs = V_r(1,:);
           new_ys = V_r(2,:);
        end
        
        function [] = predictionScoresOverTime(pf)
        % PREDICTIONSCORESOVERTIME plots prediction scores over time
        % Compares how well the AnglePlay kernel explains variance versus 
        % the maximum amount of explainable variance for a variety of 
        % offsets and winsizes.
        %
        % INPUT
        %   pf = p2m data structure
        %
        % OUTPUT
        %   (figure)
            
            winsizes = [30,50,100];
            offsets = 0:20:250;
            
            explainableVariance = zeros(length(winsizes),length(offsets));
            explainedVariance = zeros(length(winsizes),length(offsets));
            scores = zeros(length(winsizes),length(offsets));
            
            for i=1:length(winsizes)
                parfor j=1:length(offsets)
                    AP = AnglePlay(pf,offsets(j),winsizes(i));
                    explainableVariance(i,j) = AP.explainableVariance()^2;
                    explainedVariance(i,j) = AP.findKernelCorrelation(0.9,'bootstrap',1)^2;
                    scores(i,j) = explainedVariance(i,j) / explainableVariance(i,j);
                end
            end
            
            figure();
            subplot(2,2,1:2);
            plot(offsets,scores','-o');
            title([PFUtil.experName(pf) ' Kernel Prediction Scores']);
            ylabel('prediction score (explained / explainable variance)');
            xlabel('offset time (ms)');
            ylim([0 max(1,max(scores(:)))]);
            xlim([0 max(offsets)]);
            legend(cellstr(num2str(winsizes','%d ms winsize')));
            subplot(2,2,3);
            plot(offsets,explainableVariance','-o');
            title('Explainable Variance');
            ylabel('explainable variance (r^2)');
            xlabel('offset time (ms)');
            ylim([0 1]);
            xlim([0 max(offsets)]);
            legend(cellstr(num2str(winsizes','%d ms winsize')));
            subplot(2,2,4);
            plot(offsets,explainedVariance','-o');
            title('Explained Variance');
            ylabel('explained variance (r^2)');
            xlabel('offset time (ms)');
            ylim([0 1]);
            xlim([0 max(offsets)]);
            legend(cellstr(num2str(winsizes','%d ms winsize')));        
        end
        
        function [gr_gr_prediction_score,gr_ap_prediction_score,ap_ap_prediction_score] = gratrevModelComparison(ap_pf, ap_offset, ap_winsize, gr_pf, gr_offset, gr_winsize)
            
            % GratRev Gabor Model
            gr_explainable_var = RasterUtil.explainableVariance(gr_pf,gr_offset,gr_winsize)^2;
            GRF = GratRevFitter(gr_pf,gr_offset,gr_winsize);
            gr_modelfit_results = GRF.fitComplexCellSimpleGA();
            gr_gr_r_squared = gr_modelfit_results.corr;
            gr_gr_prediction_score = gr_gr_r_squared / gr_explainable_var;
                        
            % AnglePlay Kernel
            AP = AnglePlay(ap_pf,ap_offset,ap_winsize);
            ap_explainable_var = RasterUtil.explainableVariance(ap_pf,ap_offset,ap_winsize)^2;
            ap_ap_r_squared = AP.findKernelCorrelation(0.9,'bootstrap',1)^2;
            ap_ap_prediction_score = ap_ap_r_squared / ap_explainable_var;
            
            % Use Gratrev Model on AnglePlay
            load('Variables/CX_fit_args');
            load('Variables/gratrev_stimulus_space');
            CX_fit_args.orientation = gr_modelfit_results.ori;
            CX_fit_args.gabor_params.wavelength = gr_modelfit_results.wavelength;
            CX_fit_args.gabor_params.sigma = gr_modelfit_results.sigma;
            CX = ComplexCellModel(CX_fit_args);
            [~, gr_ap_rates_observed, gr_ap_rates_predicted] = AP.findModelKernel2(CX);
            gr_ap_r_squared = corr(gr_ap_rates_observed, gr_ap_rates_predicted)^2;
            gr_ap_prediction_score = gr_ap_r_squared / ap_explainable_var;
            
            figure();
            bar([gr_gr_prediction_score,gr_ap_prediction_score,ap_ap_prediction_score]);
        end
        
    end
end
