classdef PFUtil < SuperUtil
    % PFUTIL performs useful operations on p2m files
    % This class contains analyses that might be useful for any given p2m
    % file, no matter the task.
    
    properties
             
    end
    
    properties (Constant)
        FRAME_PREFIX = cell({ ...
            '2spotmap' 'spot on'; ...
            'spotmap' 'spot on'; ...
            'wnoise' 'FLIP '; ...
            'gratrev' 'FLIP '; ...
            'waverev' 'FLIP '; ...
            'ncgratrev' 'FLIP '; ...
            'msgrating' 'frameflip'; ...
            'bubbles' 'frameflip'; ...
            'reel' 'frameflip'; ...
            'flash' 'flash on'; ...
            'frzreel', 'frameflip'; ...
            'dgrat', 'SIN'; ...
            'waverev', 'FLIP '; ...
            'waverev2', 'FLIP '; ...
            'curvplay', 'frameflip '; ...
            'gridcurv', 'FLIP '
        })
    
        GAP_PREFIX = cell({ ...
            'curvplay', 'gap_frameflip ';...
            'gridcurv', 'gap_FLIP ';
        })
    
        TRIAL_PADDING = 200

    end
    
    methods

        function PFU = PFUtil()
        end
        
    end
    
    methods(Static)
       
        %% Randomization Helpers
        
        function [result_mean, result_std, results_matrix] = bootstrap(pf, f, n_draws, varargin)
        % BOOTSTRAP runs bootstrap analysis by resampling p2m file's trials
        % This function works when the function F you want to get the
        % results distribution of returns one variable. It randomly
        % resamples the trials of the experiment with replacement N_DRAWS
        % times.
            
            results = cell(n_draws,1);
            for i=1:n_draws
                rand_pf = pf;
                rand_pf.rec = rand_pf.rec(randi(length(rand_pf.rec),length(rand_pf.rec),1));
                results{i} = f(rand_pf,varargin{:});
            end
            
            result_size = size(results{1});
            results_matrix = zeros([n_draws, result_size]);
            
            for i=1:n_draws
                if min(result_size)==1
                    results_matrix(i,:) = results{i};
                else
                    results_matrix(i,:,:) = results{i};
                end
            end
            
            result_mean = squeeze(nanmean(results_matrix,1));
            result_std = squeeze(nanstd(results_matrix,1));
        end
        
        function [result_mean, result_std,results_matrix] = jackknife(pf, f, n_draws, varargin)
        % JACKKNIFE runs jackknife analysis by resampling p2m file's trials
        % This function works when the function F you want to get the
        % results distribution of returns one variable. It splits up the
        % trials into N_DRAWS chunks.        
            
            if strcmp(PFUtil.taskName(pf),'curvplay')==1
                pf = PFUtil.copyParam(pf, 'IMAGE_INFO');
            end
            
            N_RECS = length(pf.rec);
            
            shuffled_recs = randperm(N_RECS);
            recs_per_draw = floor(N_RECS / n_draws);
            
            results = cell(n_draws,1);
            max_rec_i = 0;
            for i=1:n_draws
                rand_pf = pf;
                if i~=n_draws
                    shuffled_recs_i = (max_rec_i+1):(max_rec_i + recs_per_draw);
                else
                    shuffled_recs_i = (max_rec_i+1):length(shuffled_recs);
                end
                excluded_recs = shuffled_recs(shuffled_recs_i);
                remaining_recs = setdiff(1:N_RECS, excluded_recs);
                rand_pf.rec = rand_pf.rec(remaining_recs);
                results{i} = f(rand_pf,varargin{:});
                max_rec_i = max(shuffled_recs_i);
            end
            
            result_size = size(results{1});
            results_matrix = zeros([n_draws, result_size]);
            
            for i=1:n_draws
                if min(result_size)==1
                    results_matrix(i,:) = results{i};
                else
                    results_matrix(i,:,:) = results{i};
                end
            end
            
            result_mean = squeeze(nanmean(results_matrix,1));
            result_std = squeeze(nanstd(results_matrix,1));
        end
        
        %% Low-Level P2M Helpers
        
        function name = experName(pf)
        % EXPERNAME returns name of experiment
        % e.g. romeo0284.curvplay.004
            [~,exper_name,exper_no] = fileparts(pf.src);
            name = [exper_name exper_no];
        end
        
        function taskname = taskName(pf)
        % TASKNAME returns name of experimental task
        % e.g. curvplay
              taskname = pf.rec(1).params.X_version_info{2};
              if strcmp(taskname, 'none')
                taskname = pf.rec(1).params.X_version_info{1};
              end
              [~, r] = strtok(taskname);
              [taskname, r] = strtok(r);
              taskname = strtok(taskname, '.');
        end
        
        function prefix = framePrefix(pf)
            task = PFUtil.taskName(pf);
            is = strcmp(PFUtil.FRAME_PREFIX(:,1),task);
            if sum(is) > 0
                prefix = PFUtil.FRAME_PREFIX{is,2};
            else
                prefix = 0;
            end
        end
        
        function prefix = gapPrefix(pf)
            task = PFUtil.taskName(pf);
            is = strcmp(PFUtil.GAP_PREFIX(:,1),task);
            if sum(is) > 0
                prefix = PFUtil.GAP_PREFIX{is,2};
            else
                prefix = 0;
            end
        end
        
        function [start,stop] = trialStartStop(pf, trial_i)
            
                [~,start] = PFUtil.findEvents(pf,trial_i,'stimulate_start');
                [~,stop] = PFUtil.findEvents(pf,trial_i,'stimulate_end');
                
                % Try another starter cue
                if isempty(start)
                    [~,start] = PFUtil.findEvents(pf,trial_i,'eye_start');
                    [~,stop] = PFUtil.findEvents(pf,trial_i,'eye_stop');
                end
                
                % If that also fails, then bust out this hack:
                % use the start and stop spike times with some offsets
                TRIAL_OFFSET = 2000;
                if isempty(start)
                    warning('HACK ALERT: stimulate_start and stimulate_end events missing, hacking with spike times');
                    start = min(pf.rec(trial_i).spike_times)-TRIAL_OFFSET;
                    stop = max(pf.rec(trial_i).spike_times)+TRIAL_OFFSET;
                end
        end
                
        function n = spikeCount(pf)
        % SPIKECOUNT counts all spikes from an experiment
            n = 0;
            for i=1:length(pf.rec)
                n = n + length(pf.rec(i).spike_times);
            end
        end
        
        function new_pf = copyParam(pf, key)
        % COPYPARAM copies the value of param KEY stored in the first PF rec
        % to all subsequent recs. This is required when randomly drawing
        % trials from the p2m file (e.g. copying IMAGE_INFO params when
        % jackknifing angleplay data)
        
            new_pf = pf;
            for i=2:length(new_pf.rec)
                new_pf.rec(i).params.(key) = pf.rec(1).params.(key);
            end
        end

        function new_pf = removeBadTrials(pf)
        % REMOVEBADTRIALS removes all trials that do not have the status
        % code 'C' for Correct.
            new_pf = pf;

            valid_is = [];
            for i=1:length(pf.rec)
                result = pf.rec(i).result;
                if result(1)=='C'
                    valid_is = [valid_is; i];
                end
            end
            
            if length(valid_is)==0
                error('After removing bad trials, no trials are left!');
            end
            
            new_pf.rec = pf.rec(valid_is);
        end
        
        function [ix,ts] = findEvents(pf,n,evname,exact)
        %function [ix, ts] = p2mFindEvents(pf, n, evname, exact)
        %
        % Find events prefix-matching 'evname' in the nth record of pf.
        %
        % INPUT
        %   pf = p2m data strcture
        %   n = record number
        %   evname = string to search for in event table
        %   exact = use exact matching, instead of strmatch
        %
        % OUTPUT
        %   ix = indices of matching events in pf.rec(n).ev_e/ev_t
        %   ts = times (in ms!) of matching events.
        %
        %
        % <<part of pype/p2m toolbox>>
        %
        %Thu Mar 27 22:18:31 2003 mazer 

            %pf=p2mLoad(pf);

            if ~exist('exact', 'var')
              exact = 0;
            end

            if exact
                ix_logical = strcmp(evname, pf.rec(n).ev_e);
            else
                ix_logical = strncmp(evname, pf.rec(n).ev_e, length(evname));
            end
            
            ts = pf.rec(n).ev_t(ix_logical);
            ix = find(ix_logical==1);
        end

        %% Stimulus and Rate Matrix Extractors
                
        function spikes = spikeVectors(pf, varargin)

            
            p = inputParser;
            addRequired(p,'pf');
            addParameter(p,'concatenate',0);
            addParameter(p,'padding',0); % add padding to beginning and end of trials
            addParameter(p,'bin_width', 0);
            addParameter(p,'rec_i',0);
            parse(p,pf,varargin{:});
            
            N_TRIALS = length(pf.rec);
            
            spikes = cell(N_TRIALS,1);
            
            for i=1:N_TRIALS
                
                [start,stop] = PFUtil.trialStartStop(pf,i);
                
                spike_times = pf.rec(i).spike_times;

                stop = stop - start + 1;
                spike_times = spike_times - start + 1;

                spike_times = spike_times(spike_times > 1);
                spike_times = spike_times(spike_times <= stop);
                                
                vector = zeros(stop, 1);
                vector(spike_times) = 1;
                
                spikes{i} = vector;
                
                % Add padding if concatenating trials together
                if p.Results.concatenate || p.Results.padding
                    spikes{i} = [zeros(PFUtil.TRIAL_PADDING,1); spikes{i}; zeros(PFUtil.TRIAL_PADDING,1)];
                end
            end
            
            if p.Results.concatenate
                spikes = cell2mat(spikes);
                
                if p.Results.bin_width > 0
                    spikes = SuperUtil.binVector(spikes,p.Results.bin_width);
                    spikes = spikes * p.Results.bin_width; % converts rates back into spike counts
                end
            end
        end
        
        function rates = rateVectors(pf, varargin)
            
            p = inputParser;
            addRequired(p,'pf');
            addParameter(p,'concatenate',0);
            addParameter(p,'padding',0); % add padding to beginning and end of trials
            addParameter(p,'bin_width', 0);
            addParameter(p,'rec_i',0);
            parse(p,pf,varargin{:});
            
            N_TRIALS = length(pf.rec);
            
            spikes = PFUtil.spikeVectors(pf);
            rates = cell(N_TRIALS,1);
            
            for i=1:N_TRIALS
                spike_vector = spikes{i};
                spike_density = RasterUtil.spikeDensity(spike_vector);
                rates{i} = spike_density * 1000;
                
                % Add padding if concatenating trials together
                if p.Results.concatenate || p.Results.padding
                    rates{i} = [zeros(PFUtil.TRIAL_PADDING,1); rates{i}; zeros(PFUtil.TRIAL_PADDING,1)];
                end
            end
            
            if p.Results.rec_i ~= 0
            % return only a fraction of the trials
                rates = rates(p.Results.rec_i);
            end
            
            if p.Results.concatenate
                rates = cell2mat(rates);
                
                if p.Results.bin_width > 0
                    rates = SuperUtil.binVector(rates,p.Results.bin_width);
                end
            end
        end
        
        
        function stims = stimulusMatrices(pf, varargin)
            
            p = inputParser;
            addRequired(p,'pf');
            addParameter(p,'concatenate', 0);
            addParameter(p,'bin_width', 0);
            addParameter(p,'rec_i', 0);
            addParameter(p,'sparse', 0);
            addParameter(p,'padding',0);% add padding to beginning and end of trials
            addParameter(p,'feature','binary');
            parse(p,pf,varargin{:});
            
            task = PFUtil.taskName(pf);
            frame_prefix = PFUtil.framePrefix(pf);
            gap_prefix = PFUtil.gapPrefix(pf);
            
            N_TRIALS = length(pf.rec);
            
            if strcmp(p.Results.feature,'delta')
                % do nothin bc this option is a computational walk in the
                % park!
            elseif strcmp(task,'curvplay')
                if strcmp(pf.rec(1).params.imagedir,'/auto/data/stimuli/jackson/angles_bold-0.25')
                    APL = AnglePlayLoader('imagedir','/lab/stimuli/curvplay_jackson/angles_medium/');
                else
                    APL = AnglePlayLoader();
                end
                AP = AnglePlay(pf,50,50);
                if strcmp(p.Results.feature,'binary')
                    n_oris = length(AP.orientations);
                    n_sfs = length(AP.coefficients);
                    
                    if n_sfs~=11
                        warning(['code assumes there are 11 sfs, there are actually ' num2str(n_sfs)]);
                    end
                    n_sfs_final = 6;
                    
                    % n_stim_dimensions = n_sfs*n_oris*n_types
                    if n_oris==16
                        n_stim_dimensions = n_sfs_final * 8 * 2;
                    else
                        n_stim_dimensions = n_sfs_final * n_oris * 2;
                    end
                end
                % FEATURE_TYPE = 'binary'; % 'params','binary','image'
            elseif strcmp(task,'gratrev')
                GR = GratRev(pf,50,50);
                if strcmp(p.Results.feature,'image')
                    image_size = 12;
                    n_stim_dimensions = image_size^2;
                    GRL = GratRevLoader();
                elseif strcmp(p.Results.feature,'binary')
                    n_stim_dimensions = length(GR.oris) * length(GR.sfs);
                end
            elseif strcmp(task,'waverev2')
                n_stim_dimensions = 50;
                
                gabor_list = regexp(pf.rec(1).params.ALLNAMES, '\ ', 'split');
                m = arrayfun(@(x) str2double(char(x)), ...
                       regexp(pf.rec(1).params.ALLNAMES, '[\ ,]', 'split'));
                m = reshape(m, [8 length(m)/8])';

                n_gabors = size(gabor_list,2);

                % Maps gabor strings to their index in binary vector
                gabor_map = containers.Map(gabor_list,1:n_gabors);
            end
                
            stims = cell(N_TRIALS,1);
            
            % Add Padding to Beginning and End of each Trial?
            if p.Results.concatenate || p.Results.padding
                padding = PFUtil.TRIAL_PADDING;
            else
                padding = 0;
            end
            
            % Load Raster
            unique_triggers = RasterUtil.uniqueTriggers(pf);
            stimulus_features = cell(length(unique_triggers),1);
                        
            % Create stimulus_features cell which stores a feature vector
            % corresponding to each element of the unique_triggers cell.
            % The representation of the stimulus_features vectors depends
            % on the task.
            for i=1:length(stimulus_features)
                if strcmp(p.Results.feature,'delta')
                    feature_v = [1];
                    n_stim_dimensions = 1;
                    stimulus_features{i} = feature_v;
                elseif strcmp(task,'curvplay')
                    if strcmp(p.Results.feature,'binary')
%                         n_stim_dimensions = length(pf.rec(1).params.IMAGE_INFO);
%                         stim_n = sscanf(unique_triggers{i},'%*s %d',1) + 1;
%                         feature_v = zeros(n_stim_dimensions,1);
%                         feature_v(stim_n) = 1;
                        stimParams = AnglePlayUtil.trigger2stimulusParams(pf, unique_triggers{i});
                        sf_i = find(AP.coefficients==stimParams.coefficient,1,'first');
                        ori_i = find(AP.orientations==stimParams.ori,1,'first');
                        
                        if sf_i > 1
                            sf_i = ceil((sf_i-1)/2)+1;
                        end
                        
                        if n_oris==16
                            ori_i = ceil(ori_i/2);
                        end
                        % determine binary feature index based on the following
                        % formula (collapsing across stroke):
                        % (((coeff_i - 1) * n_oris) + ori_i) + ((n_stim_dimensions/2)*type)
                        %stim_n = (((sf_i - 1) * n_oris) + ori_i) + ((n_stim_dimensions/2)*stimParams.type);
                        stim_n = (((ori_i - 1) * n_sfs_final) + sf_i) + ((n_stim_dimensions/2)*stimParams.type);
                        feature_v = zeros(n_stim_dimensions,1);
                        feature_v(stim_n) = 1;
                    elseif strcmp(p.Results.feature,'image')
                        IMG_WIDTH = 16;
                        stim_n = sscanf(unique_triggers{i},'%*s %d',1);
                        image = APL.getByImageNumber(pf, stim_n, IMG_WIDTH);
                        feature_v = reshape(image,[IMG_WIDTH^2 1]);
                        n_stim_dimensions = IMG_WIDTH^2;
                    else
                        stimParams = AnglePlayUtil.trigger2stimulusParams(pf, unique_triggers{i});
                        n_stim_dimensions = 5;
                        feature_v = [sind(stimParams.angle); cosd(stimParams.angle); sind(stimParams.ori); cosd(stimParams.angle); stimParams.type];
                    end
                    stimulus_features{i} = feature_v;
                elseif strcmp(task,'gratrev') && strcmp(p.Results.feature,'image')
                    grating = GRL.getByTrigger(pf, unique_triggers{i}, sqrt(n_stim_dimensions));
                    grating = grating - 0.5; % demean, rescaling from [0,1] -> [-0.5, 0.5]
                    stimulus_features{i} = reshape(grating,[n_stim_dimensions 1]);
                elseif strcmp(task,'gratrev') % binary
                    
                    stimParams = GratRevUtil.trigger2stimulusParams(pf, unique_triggers{i});
                    sf_i = find(GR.sfs==stimParams.sf,1,'first');
                    ori_i = find(GR.oris==stimParams.ori,1,'first');
                    n_oris = length(GR.oris);
                    % determine binary feature index based on the following
                    % formula:
                    % (((sf_i - 1) * n_oris) + ori_i)
                    stim_n = (((sf_i - 1) * n_oris) + ori_i);
                    feature_v = zeros(n_stim_dimensions,1);
                    feature_v(stim_n) = 1;
                    stimulus_features{i} = feature_v;
                elseif strcmp(task,'gridcurv') 
                    trigger = unique_triggers{i};
                    trigger = trigger(5:end); % get rid of "FLIP" prefix
                    tokens = textscan(trigger,'%*f %*f %*f %*f %f','delimiter',',');
                    oris = tokens{1};
                    
                    n_stim_dimensions = 2*length(oris);
                    feature_v = zeros(n_stim_dimensions,1);
                    for j=1:length(oris)
                        feature_v((2*(j-1)+1):(2*j)) = [sind(oris(j)); cosd(oris(j))];
                    end
                    stimulus_features{i} = feature_v;
                elseif strcmp(task,'waverev2')
                    trigger = unique_triggers{i};
                    trigger = trigger(6:end);
                    tokens = regexp(trigger,'\ ','split'); %split on spaces

                    % Loop over tokens, starting on index 2, skipping 'FLIP' token
                    wavelet_v = zeros(944,1);
                    for t=1:length(tokens)
                        gabor_i = gabor_map(tokens{t}); % look up index in dictionary
                        wavelet_v(gabor_i) = 1;
                    end
                    feature_v = wavelet_v(1:n_stim_dimensions);
                    stimulus_features{i} = feature_v;
                else
                    n_stim_dimensions = length(stimulus_features);
                    feature_v = zeros(n_stim_dimensions,1);
                    feature_v(i) = 1;
                    stimulus_features{i} = feature_v;
                end
            end
            
            for i=1:N_TRIALS
                [start,stop]=PFUtil.trialStartStop(pf,i);
                stop = stop - start + 1;
                duration = stop + 2*padding; % adds padding if going to concatenate later
                
                [frames_ix,frames] = PFUtil.findEvents(pf,i,frame_prefix);
                frames = frames - start + padding + 1;
                
                if gap_prefix ~= 0
                    [~,gaps] = p2mFindEvents(pf, i, gap_prefix);
                    gaps = gaps - start + padding + 1;
                    ends_with_gap = (length(frames) == length(gaps));
                end
                
                matrix = zeros(n_stim_dimensions,duration);
                for j=1:length(frames)
                    trigger = pf.rec(i).ev_e{frames_ix(j)};
                    trigger_i = strcmp(unique_triggers, trigger);
                    features = stimulus_features{trigger_i};
                    t_nextflip = pf.rec(i).ev_t(frames_ix(j)+1) - start + padding + 1;
                    
                    %stim_n = sscanf(,'%*s %d',1) + 1;
                    if strcmp(p.Results.feature,'delta')
                        feature_i = frames(j);
                    elseif gap_prefix ~= 0
                        if ends_with_gap || j~=length(frames)
                            feature_i = frames(j):gaps(j);
                        elseif j==length(frames)
                            feature_i = frames(j):t_nextflip;
                        end
                    else
                        if j~=length(frames)
                            feature_i = frames(j):frames(j+1);
                        elseif j==length(frames)
                            % @todo streamline if/else logic so that this
                            % fencepost case applies to both expers with
                            % and without gaps
                            feature_i = frames(j):t_nextflip;
                        else
                            error('if/else skipping cases');
                        end
                    end
                    matrix(:,feature_i) = repmat(features,[1 length(feature_i)]);
                end
                stims{i} = matrix;
            end
            
            if p.Results.rec_i ~= 0
            % return only a fraction of the trials
                stims = stims(p.Results.rec_i);
            end
                        
            if p.Results.concatenate
                
                stims = cell2mat(stims');
                
                if p.Results.bin_width > 0
                    stims = SuperUtil.binMatrix(stims,p.Results.bin_width,2);
                    
                    if strcmp(p.Results.feature,'delta')==1
                        stims(stims>0) = 1;
                        stims = stims';
                    end
                    % threshold results
                    %stims(stims < 0.5) = 0;
                    %stims(stims >= 0.5) = 1;
                end
            end
            
            if p.Results.sparse
                stims = sparse(double(stims));
            end
        end
        
        function rates = ratesCell2Vector(rates_cell,bin_width)
            rates = cell2mat(rates_cell);
            rates = SuperUtil.binVector(rates,bin_width);
        end
        
        function stims = stimulusCell2Matrix(stim_cell,bin_width)
            stims = cell2mat(stim_cell');
            stims = SuperUtil.binMatrix(stims,bin_width,2);
            if size(stims,1) > size(stims,2)
                stims = stims';
            end
        end
        
        function spikes = concatenateSpikes(pf)
        % Deprecated. Use spikeVectors() instead.
            
            n_spikes = PFUtil.spikeCount(pf);
            spikes = zeros(n_spikes,1);
            max_spikes_i = 0;

            for i=1:length(pf.rec)
                trial_spikes = pf.rec(i).spike_times;
                
                if isempty(trial_spikes)
                    continue;
                end
                
                spikes_i = (max_spikes_i + 1):(max_spikes_i+length(trial_spikes));
                spikes(spikes_i) = trial_spikes + max(spikes) + 1000;
                max_spikes_i = max(spikes_i);
            end
        end

        
        %% Stimulus Temporal Analysis
        
        function [f, mean_fdur, mean_gapdur] = stimulusCarrierFrequency(pf)
            params = pf.rec(1).params;
            taskname = PFUtil.taskName(pf);
            found_empirically = 0;
            if strcmp(taskname,'curvplay') || strcmp(taskname,'gridcurv')
                fdur = cell2mat(params.fdur);
                gapdur = cell2mat(params.gapdur);
                mean_fdur = mean(fdur);
                mean_gapdur = mean(gapdur);
            elseif strcmp(taskname,'msgrating')
                mean_fdur = params.target_hold;
                mean_gapdur = params.gap_hold;
			elseif isfield(params,'dur')
                mean_fdur = abs(params.dur);
                mean_gapdur = 0;
            else
                % if can't find via params, get an empirical estimate
                raster = prast(pf);
                isis = diff(raster.trialtimes);
                mean_fdur = median(isis);
               % warning('cannot find frame duration, using default');
               % mean_fdur = 200;
                mean_gapdur = 0;
                found_empirically = 1;
            end
            
            % if can't find via params, get an empirical estimate
            raster = prast(pf);
            isis = [];
            for i=1:length(pf.rec)
                flip_ts = raster.trialtimes(raster.trial==i);
                isis = [isis; diff(flip_ts)];
            end

            mean_fdur = mean(isis);
            % warning('cannot find frame duration, using default');
            % mean_fdur = 200;
            mean_gapdur = 0;
            found_empirically = 1;
            
            period = mean_fdur + mean_gapdur;
            f = 1 / (period/1000);
            if strcmp(taskname,'curvplay') && nargout < 1
                fprintf('fdur=[%d, %d] gapdur=[%d, %d] period=%.0fms freq=%.1fhz\n',...
                    fdur(1), fdur(2), gapdur(1), gapdur(2), period, f);
            end
        end
        
        function plotCarrierFrequency(pf, offset)
            
            if nargin < 2
                offset = 70;
            end
            offset = -offset;
            
            % PSTH
            psth = phist(pf,'pre',-1000,'post',1000);
            ts = psth(:,1);
            rates = psth(:,2);
            rates_smooth = smooth(rates,20);
            mean_psth_rate = mean(rates);
            mean_rate = SuperUtil.meanRate(pf);
            
            % Carrier Wave
            f = PFUtil.stimulusCarrierFrequency(pf);
            period = (1/f)*1000; %in ms
            angular_freq = (2*pi) / period;
            phase_px = offset; % hand-tuned, in pixels, corresponds to latency
            phase = (2 * phase_px * pi)/period;
            amplitude = mean_psth_rate / 3;%mean(psth(:,2))/2;
            y_offset = 0; %amplitude;
            carrier = amplitude*sin(angular_freq*ts + phase) + y_offset;
            
            eshade(ts,rates_smooth,psth(:,3));
            hold on;
            plot(ts,rates_smooth,'-r');
           % hold on;
            plot(ts,carrier,'-b');
            line(xlim,[0 0],'Color','k');
            set(gca,'XGrid','on');
            title('PSTH + Stimulus Carrier Frequency');
            ylabel('firing rate (s/s)');
            xlabel('time (ms)');
        end
        
        function [xs,autocorrs] = stimulusAutocorrelation(pf)
            
            vectors = PFUtil.stimulusVectors(pf);
            
            max_vector_length = 1;
            for i=1:size(vectors,1)
                ac_length = 2*length(vectors{i}) - 1;
                if ac_length > max_vector_length
                    max_vector_length = ac_length;
                end
            end
            acc = zeros(max_vector_length,1);
            
            for i=1:size(vectors,1)
                trial_autocorr = xcorr(vectors{i});
                trial_autocorr = trial_autocorr / max(trial_autocorr);
                
                l = length(trial_autocorr);
                leftover = max_vector_length - l;
                padding = floor(leftover / 2);
                acc((1 + padding):(padding + l)) = acc((1 + padding):(padding + l)) + trial_autocorr;
            end
            
            xs = -floor(max_vector_length/2):floor(max_vector_length/2);
            autocorrs = acc / size(vectors,1);
            
            figure();
            plot(xs,autocorrs);

        end
        
        function stims = stimulusVectors(pf, varargin)
        % DEPRECATED
        % Use stimulusMatrices instead. This is still used with
        % stimulusAutocorrelation().
            
            p = inputParser;
            addRequired(p,'pf');
            addParameter(p,'binary', 1);
            parse(p,pf,varargin{:});

            N_TRIALS = length(pf.rec);
            
            stims = cell(N_TRIALS,1);
            
            for i=1:N_TRIALS
                [start,stop]=PFUtil.trialStartStop(pf,i);
                [frames_ix,frames] = p2mFindEvents(pf,i,'frameflip');
                [~,gaps] = p2mFindEvents(pf,i,'gap_frameflip');
                stop = stop - start + 1;
                frames = frames - start + 1;
                gaps = gaps - start + 1;
                start = 1;
                
                ends_with_gap = (length(frames) == length(gaps));
                
                vector = zeros(stop, 1);
                vector_length_init = length(vector);
                for j=1:length(frames)
                    if ends_with_gap || j~=length(frames)
                        vector_is = frames(j):gaps(j);
                    elseif j==length(frames)
                        vector_is = frames(j):stop;
                    end
                    
                    if p.Results.binary==1
                        vector(vector_is) = 1;
                    else
                        stim_n = sscanf(pf.rec(i).ev_e{frames_ix(j)},'%*s %d',1) + 1;
                        vector(vector_is) = stim_n;
                    end
                end
                vector_length_final = length(vector);
                if vector_length_init ~= vector_length_final
                    warning('not equal!');
                end
                
                stims{i} = [zeros(PFUtil.TRIAL_PADDING,1), vector, zeros(PFUtil.TRIAL_PADDING,1)];
            end
        end

        
    end
end