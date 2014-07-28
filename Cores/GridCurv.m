classdef GridCurv < handle
    %GRIDCURV wrapper for gridcurv trial p2m files
    %   Contains utility methods to process p2m files and analyze data.
    %
    %   USAGE
    %   pf = dbfind('romeo*277*gridcurv*6');
    %   GC = GridCurv(pf, 70);
    %
    %   NOTE ON NAMING CONVENTIONS
    %   - Methods that modify the state of the object are prefixed with 'set'
    %   - Methods that do some processing on the object's properties and then 
    %        return results are prefixed with 'find'
    %   - Some methods don't fit these conventions (bootstrap, tuningStats)
    %
    %   TO-DO LIST
    %       - setFlips should not modify p2m object directly
    %       ? orilist should never be hard-coded in GridCurv methods, and
    %         instead should be extracted from the p2m object and defined
    %         as a property of the class
    %       - adopt consistent style: camelCase or under_score
    %       - add options to findOrientationMeans so that it examines
    %         certain ranges of trials at a time, or ignores frames
    %         with/without a given orientation at a position (for conditional
    %         analyses)
    %       ? find more intelligent way to set orilist
    %       - matrix sizes should never be dependent on nori
    %       - add ori @ position conditional filtering to generateStats
    %       - rename setFlips() to something more general like p2mProcess()
    %         since it also grabs grid_centers
    %       ? implement STC analysis
    %       - scale eigen vector plots by eigenvalue
    %       - move plot_gridcurv() functionality here
    
     
    properties
        
        % latency time, in ms
        offset
        
        % size of bin
        winsize
        
        % p2m object
        pf
        
        % params from p2m file
        params
        
        % shape params derived from p2m file
        shape_params
        
        % axis of grid [xmin xmax ymin ymax]
        axis
        
        % stores grid in x,y coordinate form
        grid_xy
        
        % stores grid in Axes OuterPosition form
        grid_position
        
        % list of orientations
        orientations
        
        % probability of each orientation at each position
        orientation_means
        
        % distribution of durations
        durations_dist
        
        % stimulus average
        stimulus_avg
        
        % one row per ms
        stimuli
        
        % spike-triggered stimulus ensemble
        sts

    end
    
    methods
        
        function GC = GridCurv(pf, offset, winsize, varargin)
        % GRIDCURV constructor for GridCurv class
        %   PARAM pf            (p2m object) data from gridcurv trial
        %   PARAM offset        (int) start of bin / cell's latency, in ms
        %   PARAM winsize       (int) size of bin, in ms
        %   RETURN GC           (GridCurv object)
        

            p = inputParser;
            addRequired(p,'pf');
            addRequired(p,'offset');
            addRequired(p,'winsize');
            try
              parse(p,pf,offset,winsize,varargin{:});
            catch
              error('usage: GridCurv(pf, offset, winsize, ...opts..');
            end

            GC.offset = offset;
            GC.winsize = winsize;
            GC.pf = getpf(pf);
            GC.grid_xy = [];
            
            %% Extract relevant parameters from p2m File
            GC.params = GC.pf.rec(1).params;
            
            % Calculate shape params from ratios
            GC.shape_params = {};
            GC.shape_params.rf_size = GC.params.rf_size * GC.params.scale;
            GC.shape_params.spacing_radius = (GC.params.rf_size / 2) * GC.params.spacing_ratio;
            GC.shape_params.gabor_sigma = 0.5 * GC.shape_params.spacing_radius * GC.params.sigma_ratio;
                        
            % Set orientation list
            if isfield(GC.params,'orilist') && ~isempty(GC.params.orilist)
                GC.orientations = GC.params.orilist;
            else
                if GC.params.nori==8
                    GC.orientations = [0,22,45,68,90,112,135,158];
                else
                    GC.orientations = round(0:(180/GC.params.nori):179);
                end
            end
            
            %% Process p2m File
            GC.setFlips();
            GC.setStimuli();
            GC.sts = GC.findSTS(offset, winsize);
            
            if size(GC.sts,1) < 100
                warning(['GridCurv detected only ' num2str(size(GC.sts,1)) ' spike-triggered stimuli. This may cause some analysis code to fail. Please verify that the p2m file has trials with the result code `C Correct`.']);
            end
            [GC.orientation_means, GC.stimulus_avg] = GC.findOrientationMeans();
            
            %% Normalize plotPolars grid coordinates into axes OuterPosition format
            if size(GC.sts,2)==1
                GC.axis = [-10 10 -10 10];%min(GC.grid_xy(:,1)) max(GC.grid_xy(:,1)) min(GC.grid_xy(:,2)) max(GC.grid_xy(:,2))];
            else
                GC.axis = [min(GC.grid_xy(:,1)) max(GC.grid_xy(:,1)) min(GC.grid_xy(:,2)) max(GC.grid_xy(:,2))];
            end

            if size(GC.sts,2)==1
                GRID_WIDTH = 1;
                GRID_HEIGHT = 1;
                GRID_SCALE = 1;
            elseif size(GC.sts,2)==7
                GRID_WIDTH = 0.35;
                GRID_HEIGHT = 0.3;
                GRID_SCALE = 1.5;
            elseif size(GC.sts,2)==19
                GRID_WIDTH = 0.24;
                GRID_HEIGHT = 0.2;
                GRID_SCALE = 1.3;
            end
            
            % left, bottom, width, height
            GC.grid_position(:,1) = (GC.grid_xy(:,1)+abs(GC.axis(1)))/((abs(GC.axis(1))+abs(GC.axis(2)))*GRID_SCALE);
            GC.grid_position(:,2) = (GC.grid_xy(:,2)+abs(GC.axis(3)))/((abs(GC.axis(3))+abs(GC.axis(4)))*GRID_SCALE);
            GC.grid_position(:,3) = GRID_WIDTH;
            GC.grid_position(:,4) = GRID_HEIGHT;
            
            % if only one cell, display fullscreen
            if(size(GC.grid_position,1)==1)
                GC.grid_position(1,:) = [0,0,GRID_WIDTH,GRID_HEIGHT];
            end

        end
        
        function [] = setFlips(GC)
        % SETFLIPS computes frame flip times from GC.pf
        %   PARAM GC            (object) implicit self parameter

            % Flip event entries are encoded as a {N_CELL} row, 5 column table
            %    1: CIRCLE INDEX
            %    2: CENTER X
            %    3: CENTER Y
            %    4: PHASE
            %    5: ORI

            n_flips = 0;
            
            if GC.params.nlevels==0
                n_cells = 1;
            elseif GC.params.nlevels==1
                n_cells = 7;
            elseif GC.params.nlevels==2
                n_cells = 19;
            elseif GC.params.nlevels==3
                n_cells = 37;
            end
                

            %% Loop over each trial, converting flips
            for r=1:size(GC.pf.rec,2)

                flip_i = 1;

                % Create copy of event info which will store the vector format
                % The flips cell is a n x 4 table which stores the following info:
                %	{n, 1} = event time
                % 	{n, 2} = all grid cell entries as a single string
                % 	{n, 3} = grid cells as a 7 x 1 orientation vector
                flips = cell(1,3);

                % flips_t vector stores list of times when flips occurred, 
                % and is used to lookup stimuli by spike timestamp
                flips_t = [];
                gaps_t = [];

                %% Loop over each event in the trial and convert into binary vector
                for i=1:max(size(GC.pf.rec(r).ev_e))
                    % If the event code starts with 'FLIP', we know that event with
                    % index of i is a flip event that is to be processed
                    if strncmp(GC.pf.rec(r).ev_e{i},'FLIP',4)==1

                        % Copy time to new array
                        flips{flip_i,1} = GC.pf.rec(r).ev_t(i);
                        flips_t = [flips_t ; GC.pf.rec(r).ev_t(i)];

                        % Copy string to new array
                        flips{flip_i,2} = GC.pf.rec(r).ev_e{i};
                        
                        %%  Convert gabor string to binary
                        grid_bin = zeros(n_cells,1,'uint8');
                        cells = regexp(GC.pf.rec(r).ev_e{i},'\ ','split'); %split on spaces

                        % Loop over tokens, starting on index 2, skipping 'FLIP' token
                        for c=2:max(size(cells))

                            params = regexp(cells(c),',','split');
                            params = params{1};
                            grid_bin(str2double(params{1})+1) = str2double(params{5});
                        end
                        
                        % Grab grid centers (if not yet set)
                        if isempty(GC.grid_xy)
                            for c=2:max(size(cells))
                                params = regexp(cells(c),',','split');
                                params = params{1};
                                GC.grid_xy(str2double(params{1})+1,:) = [str2double(params{2}),str2double(params{3})];
                            end
                        end
                        
                        flips{flip_i,3} = grid_bin;

                        % increment counter
                        flip_i = flip_i + 1;
                    elseif strncmp(GC.pf.rec(r).ev_e{i},'gap',3)==1
                        gaps_t = [gaps_t ; GC.pf.rec(r).ev_t(i)];
                    end
                end

                n_frames = size(flips,1);

                GC.pf.rec(r).flips = flips;
                GC.pf.rec(r).flips_t = flips_t;
                GC.pf.rec(r).gaps_t = gaps_t;

                n_flips = n_flips + flip_i;
            end

            GC.pf.flips = 1;
            GC.pf.n_flips = n_flips;
        end
        
        function [] = setStimuli(GC)
        % SETSTIMULI extracts stimuli vectors at 1 ms resolution
        %   Sets GC.stimuli and GC.durations_dist
        %
        %
        %   PARAM GC                (object) implicit self parameter

            % If no ms_max, define to be impossibly large.
            N_CELLS = size(GC.grid_xy,1);
            ms_i = 1;
            frames_i = 1;            
            GC.stimuli = zeros(50000*size(GC.pf.rec,2),N_CELLS);
            GC.durations_dist = zeros(100*size(GC.pf.rec,2),1);
            
            % Loop over each trial
            for r=1:size(GC.pf.rec,2)

                % Only process correct/valid trials
                if strcmp(GC.pf.rec(r).result,'C Correct')==0
                    continue;
                end
                
                %% Loop over each flip
                for i=1:size(GC.pf.rec(r).flips)

                    % Plan:
                    % 1) get flip time
                    % 2) get next gap time
                    % 2b) if no next gap found, get stimulate_end time
                    %     verify that stimulate_end - flip < 500

                    stim_start = GC.pf.rec(r).flips{i,1};

                    % No More Gaps! Find stimulate_end
                    if max(GC.pf.rec(r).gaps_t > stim_start)==0
                        [~,t] = PFUtil.findEvents(GC.pf,r,'stimulate_end');
                        stim_end = t;
                    else
                        gap = find(GC.pf.rec(r).gaps_t > stim_start,1,'first');
                        stim_end = GC.pf.rec(r).gaps_t(gap);
                    end

                    stim_dur = stim_end - stim_start;
                    
                    % create stim vector row for each ms of stimulus duration
                    GC.stimuli(ms_i:(ms_i+(stim_dur-1)),:) = repmat(GC.pf.rec(r).flips{i,3}',stim_dur,1);
                    GC.durations_dist(frames_i) = stim_dur;

                    frames_i = frames_i + 1;
                    ms_i = ms_i+(stim_dur-1);
                end
            end

            % trim thetas and durations matrices
            GC.stimuli = GC.stimuli(1:ms_i,:) - GC.params.bestori;
            GC.durations_dist = GC.durations_dist(1:(frames_i-1),:);
        end
        
        function [spike_triggered_stimuli] = findSTS(GC, offset, winsize)
        % FINDSTS searches GC.pf for the spike-triggered-stimuli by examining
        % the spikes in the window (with duration WINSIZE in ms) some OFFSET
        % ms after stimulus onset. 
        %
        % 	Note: GC.pf must first be processed using setFlips()!
        %
        %   PARAM GC        (object) implicit self parameter
        %   PARAM offset    (int) lag time in ms
        %   PARAM winsize   (int) window size in ms
        %   RETURNS         (n_spikes x n_oris matrix) spike_triggered_stimuli

            spike_triggered_stimuli = zeros(1000*size(GC.pf.rec,2),size(GC.grid_xy,1),'uint8');
            spike_i = 1;

            % Loop over each trial
            for r=1:size(GC.pf.rec,2)

                % Only process correct/valid trials
                if strcmp(GC.pf.rec(r).result,'C Correct')==0
                    continue;
                end

                n_spikes = max(size(GC.pf.rec(r).spike_times));

                %% Loop over each spike
                for i=1:n_spikes
                        
                    % use offset and winsize to focus on binned data
                    t_spike = GC.pf.rec(r).spike_times(i);

                    % find indices of adjacent onset times
                    stimulus_prev_i = find(GC.pf.rec(r).flips_t < t_spike,1,'last');
                    gap_prev_i = find(GC.pf.rec(r).gaps_t < t_spike,1,'last');
                    if size(stimulus_prev_i,1)==0 || stimulus_prev_i==0
                        continue;
                    end


                    % find stimulus and gap onset times
                    t_stimulus_prev = GC.pf.rec(r).flips_t(stimulus_prev_i);

                    if ~isempty(gap_prev_i)
                        t_gap_prev = GC.pf.rec(r).gaps_t(gap_prev_i);
                    else
                        t_gap_prev = 0;
                    end                       

                    % Case 1: Spike During Gap
                    if t_gap_prev > t_stimulus_prev
                        continue;

                    % Case 2: Spike During Stimulus, Before Bin
                    elseif t_spike < (t_stimulus_prev + offset)
                        continue;

                    % Case 3: Spike During Stimulus, After Bin
                    elseif t_spike > (t_stimulus_prev + offset + winsize)
                        continue;
                    end

                    % Case 4: Spike During Stimulus and During Bin
                    % Spike added to STS below...

                    spike_triggered_stimuli(spike_i,:) = GC.pf.rec(r).flips{stimulus_prev_i,3};

                    spike_i = spike_i + 1;
                end
            end
            spike_triggered_stimuli = spike_triggered_stimuli(1:(spike_i-1),:) - GC.params.bestori;
        end
        
        function [control_stimuli] = findControlStimuli(GC)
            indices = randi(size(GC.stimuli,1),size(GC.sts,1),1);
            control_stimuli = zeros(size(GC.sts));
            control_stimuli = GC.stimuli(indices,:);
        end
        
        function [p_theta, stimulus_avg] = findOrientationMeans(GC, varargin)
        % FINDORIENTATIONMEANS finds the probability that each orientation will be
        % present at each grid position
        %
        %   This function is used to correct for any orientation
        %   anisotropy, and to normalize probability distributions.
        %
        %   PARAM GC                (object) implicit self parameter
        %   PARAM varargin.ms_max   (int) stop after ms_max milliseconds
        %   RETURNS p_theta         (n_cells x n_oris matrix) probability of each 
        %                           orientation at each position

            p = inputParser;
            addRequired(p,'GC');
            addParamValue(p,'ms_max',[]);
            parse(p,GC,varargin{:});
            
            % If no ms_max, define to be impossibly large.
            MS_MAX = p.Results.ms_max;
            N_CELLS = size(GC.grid_xy,1);
            
            % Get thetas
            if ~isempty(MS_MAX)
                thetas = GC.stimuli(1:MS_MAX,:);
            else
                thetas = GC.stimuli;
            end

            % Create p_theta average
            hist_results = zeros(N_CELLS,length(GC.orientations));
            for i=1:N_CELLS
                [counts,~] = hist(thetas(:,i),GC.orientations);
                hist_results(i,:) = counts;
            end
            p_theta = hist_results / size(thetas,1);
            
            % Create stimulus average in XY space
            thetas_xy = zeros([size(thetas,1) N_CELLS*2]);
            thetas_xy(:,1:N_CELLS) = sind(2*thetas);
            thetas_xy(:,(N_CELLS + 1):(2 * N_CELLS)) = cosd(2*thetas);
            stimulus_avg = mean(thetas_xy,1);
        end
        
        function [eigenvalues_mean,eigenvalues_std] = findControlEigenvalues(GC)
            if size(GC.sts,1) < 10000
                N_DRAWS = 40;
            else
                N_DRAWS = 20;
            end
            eigenvalues = zeros(N_DRAWS,size(GC.grid_xy,1)*2);
            for i=1:N_DRAWS
                controls = GC.findControlStimuli();
                result = GC.stcAnalysis(controls,GC.stimulus_avg);
                eigenvalues(i,:) = result.eigen_values;
            end
            eigenvalues_mean = mean(eigenvalues,1);
            eigenvalues_std = std(eigenvalues,1);
        end
        
        function [ results ] = stcAnalysis(GC, spike_stimuli, control_average)
        % STCANALYSIS Runs spike-triggered covariance analysis on data.
        %   PARAM spike_stimuli             spike-triggered stimuli
        %   PARAM control_average           avg of control stimuli, used when non-white
        %   RETURN results

            N_CELLS = size(spike_stimuli,2);
            N_SPIKES = size(spike_stimuli,1);
            N_EIGENS = N_CELLS * 2;
            
            if size(spike_stimuli,2)==N_CELLS
                spike_stimuli = double(spike_stimuli);
                spike_stimuli_orig = spike_stimuli;
                spike_stimuli_v = zeros([N_SPIKES N_CELLS*2]);
                spike_stimuli_v(:,1:N_CELLS) = sind(2*spike_stimuli);
                spike_stimuli_v(:,(N_CELLS + 1):(2 * N_CELLS)) = cosd(2*spike_stimuli);
                spike_stimuli = spike_stimuli_v;
            end
            
            biases = (1/length(GC.orientations))*ones(size(GC.orientation_means)) ./ GC.orientation_means;
            biases = ones(size(GC.orientation_means));
            
            for spike=1:spike_stimuli_orig
                for cell=1:N_CELLS
                    ori_i = find(GC.orientations==spike_stimuli_orig(spike,cell));
                    bias = biases(cell,ori_i);
                    spike_stimuli(spike,cell) = spike_stimuli(spike,cell)*bias;
                    spike_stimuli(spike,N_CELLS+cell) = spike_stimuli(spike,cell)*bias;
                end
            end
            
            %% STA
            sta = mean(spike_stimuli,1);

            % If control average supplied, subtract it from STA
            if ~isempty(control_average)
            %    sta = sta - control_average;
            end

            %% STC

            % Normalize spike_stimuli vector
            spike_stimuli_normalized = zeros(size(spike_stimuli),'single');
%             % Original Method (from CX analysis)
%             row_wise_mean = mean(spike_stimuli,2);
%             row_wise_std = std(spike_stimuli');

%             sta_n = (sta - mean(sta(:)))/std(sta(:));
%             for i=1:N_SPIKES
%                 sts_n = (spike_stimuli(i,:) - row_wise_mean(i)) / row_wise_std(i);
%                 spike_stimuli_normalized(i,:) = sts_n - sta_n;
%             end
            
            % Method #2: No normalization, still subtracting STA
            for i=1:N_SPIKES
                spike_stimuli_normalized(i,:) = spike_stimuli(i,:) - sta;
            end

%             % Method #3: No normalization, no STA subtraction
%             spike_stimuli_normalized = spike_stimuli;
            

            % Calculate covariance matrix
            st_cov = cov(spike_stimuli_normalized);

            % Principal Component Analysis
            % grab all eigen vectors using function for dense pca
            [eigen_vectors,eigen_values_diag] = eig(st_cov);
            eigen_values = diag(eigen_values_diag);

            % everything is sorted in ascending order, so we need to resort it
            % to match up with the sorting of the eigs function used below
            [eigen_values,I] = sort(eigen_values,'descend');
            eigen_vectors = eigen_vectors(:, I);

            %% Package results
            results = {};
            results.cov = st_cov;
            results.sta = sta;
            results.eigen_vectors = eigen_vectors;
            results.eigen_values = eigen_values;
        end
        
        function [stats] = generateStats(GC, varargin)
        % GENERATESTATS bootstrap analysis of gridcurv spike-triggered stimuli
        %   Bootstrap analysis run and then raw results and 1st order stats stored
        %   in the STATS object, which is then passed on to plotting functions.
        %
        %   To analyze a fraction of the dataset: provide the optional 
        %   VARARGIN.FRACTION parameter. For example, if FRACTION = 0.25 then 
        %   GENERATESTATS examines the first quarter of the dataset.
        %
        %   To manually define the sts: provide VARARGIN.STS
        %   
        %   PARAM GC                (object) implicit self parameter
        %   PARAM varagin.fraction  (double) fraction of dataset to examine
        %   PARAM varagin.sts       (n_spike x n_grid_cell matrix) portion of
        %                           sts to be analyzed
        
            p = inputParser;
            addRequired(p,'GC');

            addParamValue(p,'fraction',1);
            addParamValue(p,'sts',100);

            parse(p,GC,varargin{:});
            
            % Analyze fraction of dataset
            if p.Results.fraction == 1
                orientation_means_slice = GC.orientation_means;
                sts_slice = GC.sts;
            else
                fraction = p.Results.fraction;
                duration = sum(GC.durations_dist);
                max_ms = round(duration * fraction);
                orientation_means_slice = GC.findOrientationMeans('ms_max',max_ms);
                sts_slice = GC.sts(1:round(size(GC.sts,1)*fraction),:);
            end
            
            % Manually define STS
            % note: unless fraction also provided, uses orientation means
            % from the entire data set

            if p.Results.sts ~= 100
                sts_slice = p.Results.sts;
            end

            N_DRAWS = 250;
            N_ORI = length(GC.orientations);
            N_CELLS = size(sts_slice,2);
            stats = {};

            % Calculate orientation bias. If all orientations appeared
            % equally as often, then this will be a matrix of ones
            stats.orientation_bias = (1/N_ORI)*ones(size(orientation_means_slice)) ./ orientation_means_slice;

            % Run bootstrap, correcting for orientation biases
            bootfun = @(sts_slice) GridCurv.findTuningStats(sts_slice, GC.orientations, stats.orientation_bias);
            [stats.raw,~] = bootstrp(N_DRAWS,bootfun,sts_slice);

            % Run booststrap again, without correcting for orientation biases
            bootfun2 = @(sts_slice) GridCurv.findTuningStats(sts_slice, GC.orientations, ones(size(stats.orientation_bias)));
            [stats.raw_uncorrected,~] = bootstrp(N_DRAWS,bootfun2,sts_slice);
            
            % Calculate mean vectors from p_thetas
            vectors = zeros(N_DRAWS, N_CELLS, 4);
            for i=1:N_DRAWS
                draw_p_thetas = reshape(stats.raw(i,:)',N_CELLS, N_ORI);
                draw_vectors = GC.pThetasToMeanVectors(draw_p_thetas);
                vectors(i,:,:) = draw_vectors;
            end
            stats.vectors = squeeze(mean(vectors,1));
            stats.vectors_se = squeeze(std(vectors,0,1));

            % Compute bootstrap 1st order stats
            stats.p_theta = reshape(mean(stats.raw)',N_CELLS, N_ORI);
            stats.p_theta_uncorrected = reshape(mean(stats.raw_uncorrected)',N_CELLS, N_ORI);
            stats.se = reshape(std(stats.raw)',N_CELLS, N_ORI);
            stats.p_theta_min = reshape(prctile(stats.raw,1)',N_CELLS, N_ORI);
            stats.p_theta_max = reshape(prctile(stats.raw,99)',N_CELLS, N_ORI);
        end
        
        function [stats, neighborhood_stats] = generateConditionalStats(GC, subunit_i, ori_i)

            if ischar(subunit_i)
                if strcmp(subunit_i,'center')==1
                    % find center
                    [~,subunit_i] = max(sum(GC.grid_xy(:,:)==0,2));
                else
                    error('bad subunit_i: acceptable values are integers or string ''center''');
                end
            end
                

            ind = find(GC.sts(:,subunit_i)==GC.orientations(ori_i));
            sts = GC.sts(ind,:);

            stats = GC.generateStats('sts',sts);

            % Max +/- 1 orientation step
            neighbor_left = ori_i - 1;
            neighbor_right = ori_i + 1;
            if neighbor_left==0
                neighbor_left = length(GC.orientations);
            end
            if neighbor_right > length(GC.orientations)
                neighbor_right = 1;
            end

            ind1 = find(GC.sts(:,subunit_i)==GC.orientations(neighbor_left));
            ind2 = find(GC.sts(:,subunit_i)==GC.orientations(ori_i));
            ind3 = find(GC.sts(:,subunit_i)==GC.orientations(neighbor_right));
            ind = unique([ind1; ind2; ind3]);
            neighborhood_sts = GC.sts(ind,:);

            neighborhood_stats = GC.generateStats('sts',neighborhood_sts);

        end
        
        function [v] = pThetasToMeanVectors(GC, p_thetas)
        % PTHETASTOMEANVECTORS generates mean vectors
        %   PARAM p_thetas
        %   RETURN v            vector of angle, magnitude, x & y components
        
            N_CELLS = size(GC.sts,2);
            v = zeros(N_CELLS,4); % stores angle and magnitude vectors
            for cell=1:N_CELLS
                xs = 0;
                ys = 0;
                for theta=1:length(GC.orientations)
                    r = p_thetas(cell,theta);
                    xs = xs + (r * cosd(2*GC.orientations(theta)));
                    ys = ys + (r * sind(2*GC.orientations(theta)));
                end
                
                % convert back [0:180]
                angle = atan2d(ys,xs) / 2;
                magnitude = hypot(xs,ys);
                new_x = cosd(angle) * magnitude;
                new_y = sind(angle) * magnitude;
                v(cell,:) = [angle; magnitude; new_x; new_y];
            end
        end
        
        function [uiParent, stats] = parsePlots(GC,varargin)
        %PARSEPLOTS gets relevant params for all plot functions
        %   PARAM GC                (object) implicit self parameter
        %   PARAM varargin          (cell) contains parameter key/val pairs
        %   PARAM varargin.stats    (object) bootstrap results
        %   PARAM varargin.Parent   (handle to UI panel / figure)
        %   RETURN uiParent         (handle to UI panel / figure)
        %   RETURN stats            (object) bootstrap results
            p = inputParser;
            addRequired(p,'GC');
            addParamValue(p,'stats',0);
            addParamValue(p,'Parent',0);
            parse(p,GC,varargin{:});
            
            % set parent
            if p.Results.Parent ~= 0
                uiParent = p.Results.Parent;
            else
                uiParent = gca;
            end
            
            % set stats
            if isstruct(p.Results.stats)
                stats = p.Results.stats;
            else
                stats = GC.generateStats();
            end
        end
        
        function [] = plotBars(GC,varargin)
        % PLOTBARS plots bars
        %   PARAM GC
        %   PARAM varargin.stats        (object) bootstrap results
        %   PARAM varargin.Parent       (handle to UI panel / figure)
            
            [uiParent, stats] = GC.parsePlots(varargin{:});
            
            P_NULL = 1 / length(GC.orientations);
            
            if size(GC.sts,2)==1
                LENGTH = (GC.axis(2) - GC.axis(1));
                WIDTH = 7;
            elseif size(GC.sts,2)==7
                LENGTH = (GC.axis(2) - GC.axis(1))/3;
                WIDTH = 5;
            elseif size(GC.sts,2)==19
                LENGTH = (GC.axis(2) - GC.axis(1)) / 5.5;
                WIDTH = 4;
            elseif size(GC.sts,2)==37
                LENGTH = (GC.axis(2) - GC.axis(1)) / 8;
                WIDTH = 3;
            end

            a = get(uiParent);
            if strcmp(a.Type,'axes')
                axes(uiParent);
            else
                axes('Parent',uiParent);
                colorbar();
            end
            
            colors = hot(100);

            for i=1:size(stats.p_theta,1)

                % each row of the angles vector stores angle and weight for each 
                % bar that will be plotted at this location. for 7-vector inputs, 
                % angles has only one row. for inputs with histogram counts for
                % each angle, then angles will have one row per angle.
                angles = zeros(1,2); % angles(i,1)==ori, angles(i,2)==weight [0:1]

                % many orientations per position (i.e. histogram-count inputs)
                for j=1:length(GC.orientations)
                    angles(j,1) = GC.orientations(j);
                    angles(j,2) = stats.p_theta(i,j);
                end

                % sort so that the bars with highest weight are always on top
                weights = angles(:,2);
                [~, I] = sort(weights);
                angles = angles(I,:);

                for a=1:size(angles,1)

                    angle = angles(a,1)+90; % correct for weird gabor ori conventions

                    dx = (LENGTH/2)*cosd(angle);
                    dy = (LENGTH/2)*sind(angle);

                    xs = [(GC.grid_xy(i,1)-dx) (GC.grid_xy(i,1)+dx)];
                    ys = [(GC.grid_xy(i,2)-dy) (GC.grid_xy(i,2)+dy)];

                    weight = round(3*angles(a,2)*100);
                    if weight >= 100
                        c = colors(100,:);
                    elseif weight <= 0
                        c = colors(1,:);
                    else
                        c = colors(weight,:);
                    end

                    line(xs,ys,'LineWidth',WIDTH,'Color',c);
                    hold on;
                end
            end
            set(gca,'Color',colors(round(3*P_NULL*100),:));
            axis(1.5*GC.axis);
            set(gca,'XTick',[]); set(gca,'YTick',[]);
            axis square;
            colormap(colors);
            caxis([0 0.333]);
            
        end
        
        function [] = plotTuningCurves(GC, varargin)
        % PLOTTUNINGCURVES plots panel of tuning curves for each grid cell
        %   PARAM GC
        %   PARAM varargin.stats        (object) bootstrap results
        %   PARAM varargin.Parent       (handle to UI panel / figure)
            
            [uiParent, stats] = GC.parsePlots(varargin{:});
            
            p_null = 1 / length(GC.orientations);
            
            % Rotate so that 0 degrees is at the center
            oris = GC.orientations;
            circshift_amount = round(length(GC.orientations)/2);
            oris = circshift(oris,[0 circshift_amount]);
            oris(1:circshift_amount) = oris(1:circshift_amount) - 180;
            stats.p_theta = circshift(stats.p_theta,[0 circshift_amount]);
            stats.se = circshift(stats.se,[0 circshift_amount]);
            stats.p_theta_min = circshift(stats.p_theta_min,[0 circshift_amount]);
            stats.p_theta_max = circshift(stats.p_theta_max,[0 circshift_amount]);
   
            % copy max value to min, to make plot symmetrical
            oris = [oris, min(oris)+180];
            stats.p_theta(:,9) = stats.p_theta(:,1);
            stats.se(:,9) = stats.se(:,1);
            stats.p_theta_min(:,9) = stats.p_theta_min(:,1);
            stats.p_theta_max(:,9) = stats.p_theta_max(:,1);
            
            y_max = max(2*p_null,max(stats.p_theta_max(:)*1.1));

            % loop over each grid cell
            for i=1:size(stats.p_theta,1)
                axes('Parent',uiParent,'OuterPosition',GC.grid_position(i,:));

                errorbar(oris,stats.p_theta(i,:),stats.se(i,:),'r');
                hold on;
                plot(oris,stats.p_theta_min(i,:),'-','Color',[0.7 0.7 0.7]);
                plot(oris,stats.p_theta_max(i,:),'-','Color',[0.7 0.7 0.7]);
                plot(oris,p_null*ones(size(oris)),'-b');

                set(gca,'XTick',oris(1:2:9));
                axis([min(oris) max(oris) 0 y_max]);
                ylabel('p(\theta | spike)');
            end
        end
        
        function [] = plotPolars(GC, varargin)
        % PLOTTUNINGCURVES panel of tuning polar plots for each grid cell
        %   PARAM GC
        %   PARAM varargin.stats        (object) bootstrap results
        %   PARAM varargin.Parent       (handle to UI panel / figure)
            
           [uiParent, stats] = GC.parsePlots(varargin{:});
            
           for i=1:size(stats.p_theta,1)

                % each row of the angles vector stores angle and weight for each 
                % bar that will be plotted at this location. for 7-vector inputs, 
                % angles has only one row. for inputs with histogram counts for
                % each angle, then angles will have one row per angle.
                angles = zeros(1,4); % angles(i,1)==ori, angles(i,2)==weight [0:1]

                % many orientations per position (i.e. histogram-count inputs)
                for j=1:length(GC.orientations)
                    angles(j,1) = GC.orientations(j);
                    angles(j,2) = stats.p_theta(i,j);
                    angles(j,3) = stats.p_theta_min(i,j);
                    angles(j,4) = stats.p_theta_max(i,j);
                end

                % expand from 0-158 degrees to 0-360 degrees
                angles(9:16,1) = GC.orientations + 180;
                angles(9:16,2:4) = angles(1:8,2:4);
                angles(17,:) = angles(1,:);

                a = axes('Parent',uiParent,'OuterPosition',GC.grid_position(i,:));        

                % hack: create phantom plot with 0-1 radius, and then hide it
                max_r = 1;
                t = 0 : .01 : 2 * pi;
                P = polar(t, max_r * ones(size(t)));
                set(P, 'Visible', 'off');
                set(findall(gca, 'String', '  1'),'String', '  0.25'); 
                set(findall(gca, 'String', '  0.5'),'String', '  0.125'); 
                hold on;

                h = polar(degtorad(angles(:,1)),angles(:,2)*4,'-b','Parent',a);
                set(h,'LineWidth',2);

                h = polar(degtorad(angles(:,1)),angles(:,3)*4,'-r','Parent',a);
                set(h,'LineWidth',1);

                h = polar(degtorad(angles(:,1)),angles(:,4)*4,'-r','Parent',a);
                set(h,'LineWidth',1);

                view(-90,90);
           end
        end
        
        function [] = plotRF(GC)
        % PLOTRF plots receptive field and grid cell sigmas
        %   Circle dimensions taken from GC.shape_params
        %   dotted circle == RF
        %   solid circles == grid cells
            
            plot_circle = @(x,y,diam,style) rectangle('Position',[x - diam/2, y - diam/2, diam, diam], 'Curvature',[1,1],'LineStyle',style);
            plot_circle(0,0,GC.shape_params.rf_size,'--');
            for i=1:size(GC.grid_xy,1)
                plot_circle(GC.grid_xy(i,1), GC.grid_xy(i,2), 2*GC.shape_params.gabor_sigma,'-');
            end
            title('RF and Grid Cells Positions');
            axis square;
        end
        
        function [] = plotVectorGabors(GC, v, varargin)
        % PLOTVECTORGABORS plots gabors of individual N_CELL-length vectors
            
            p = inputParser;
            addRequired(p,'GC');
            addRequired(p,'v');
            addParamValue(p,'Parent',0);
            parse(p,GC,v,varargin{:});
            
            % set parent
            if p.Results.Parent ~= 0
                uiParent = p.Results.Parent;
            else
                uiParent = gca;
            end
            
            N_CELLS = size(GC.sts,2);
            
            weights = ones(N_CELLS,1);

            % convert vector with 2*N_CELLS x and y components into one with
            % N_CELLS angular components
            if max(size(v))==2*N_CELLS
                weights = sqrt(v(1:N_CELLS).^2 + v((N_CELLS + 1):2*N_CELLS).^2);
                %weights = weights * (1/max(weights(:))); % scale to range [0:1]
                v = atan2d(v(1:N_CELLS),v((N_CELLS + 1):2*N_CELLS));
                v = v/2;
            end

            %% Plot as Gabors

            % 65x65 Kernel Parameters
            gabor_wavelength = 25;
            gabor_phase = pi/2;
            gabor_sigma = 10;
            gabor_aspect = 1;
            
            % Image size
            img_width = max((GC.axis(2)-GC.axis(1))*1.5,max((GC.axis(4)-GC.axis(3)*1.5),200));
            if mod(img_width,2)==1
                img_width = img_width+1;
            end
            img_size = [img_width img_width];

            % Plot subunits
            subunits = zeros(img_size);
            for i=1:size(GC.grid_xy,1)
                gb = GaborUtil.createFilter(v(i),gabor_wavelength,gabor_sigma,gabor_phase,gabor_aspect);
                gb = GaborUtil.translate(gb, img_size, GC.grid_xy(i,1), GC.grid_xy(i,2));
                subunits = subunits + weights(i)*gb;
            end
            imagesc(subunits,'Parent',uiParent);
            set(gca,'YTick',[]); set(gca,'XTick',[]); axis square;
        end
        
        function [] = plotSTC(GC)
        % @todo
        end
        
        function [] = plotSubspace(GC,indices)
        % PLOTSUBSPACE plots sts in theta-space
        %   PARAM indices            vector of three subunits to display

            N_POINTS = 4000;
            sts_i = randi(size(GC.sts,1),N_POINTS,1);
            stim_i = randi(size(GC.stimuli,1),N_POINTS,1);
            
            sts_slice = GC.sts(sts_i,indices);
            stim_slice = GC.stimuli(stim_i,indices);
            
            % add noise!
            sts_noise = randn(size(sts_slice)) .* 7;
            stim_noise = randn(size(stim_slice)) .* 7;
            sts_slice = double(sts_slice) + sts_noise;
            stim_slice = double(stim_slice) + stim_noise;
            
            
            figure();
            
            subplot(2,2,1);
            %scatter3(stim_slice(:,1),stim_slice(:,2),stim_slice(:,3),5,'MarkerFaceColor','k');
            %hold on;

            scatter3(sts_slice(:,1),sts_slice(:,2),sts_slice(:,3),7,'MarkerFaceColor','r','MarkerEdgeColor','r');
            axis([-10 170 -10 170 -10 170]); axis square;
            
            subplot(2,2,2);
            title([num2str(indices(1)) ' vs ' num2str(indices(2))]);
            scatter(sts_slice(:,1),sts_slice(:,2),5);
            xlabel(['subunit ' num2str(indices(1))]);
            ylabel(['subunit ' num2str(indices(2))]);
            axis([-10 170 -10 170]); axis square;
            
            subplot(2,2,3);
            title([num2str(indices(2)) ' vs ' num2str(indices(3))]);
            scatter(sts_slice(:,2),sts_slice(:,3),5);
            xlabel(['subunit ' num2str(indices(2))]);
            ylabel(['subunit ' num2str(indices(3))]);
            axis([-10 170 -10 170]); axis square;
            
            subplot(2,2,4);
            title([num2str(indices(1)) ' vs ' num2str(indices(3))]);
            scatter(sts_slice(:,1),sts_slice(:,3),5);  
            xlabel(['subunit ' num2str(indices(1))]);
            ylabel(['subunit ' num2str(indices(3))]);
            axis([-10 170 -10 170]); axis square;
        
        end
        
        function [] = plotSubspaces(GC)
            GC.plotSubspace([1,2,4]);
            GC.plotSubspace([4,5,7]);
            GC.plotSubspace([3,4,6]);
        end
        
        
        function [] = plotCrossValidation(GC,varargin)
        % PLOTCROSSVALIDATION plots scores as function of number of trials
        %
        % This function is useful for understanding when enough data has
        % been collected.
            
            p = inputParser;
            addRequired(p,'GC');
            addParamValue(p,'Parent',0);
            parse(p,GC,varargin{:});
            
            % set parent
            if p.Results.Parent ~= 0
                uiParent = p.Results.Parent;
            else
                uiParent = gca;
            end
            
            % Measure Ground Truth (aka full-dataset)
            full_stats = GC.generateStats();
            truth = full_stats.p_theta(:);
            
            % Generate Fractional P-Thetas
            fractions = 0.2:0.1:1;
            p_thetas = cell(length(fractions),1);
            correlations = zeros(length(fractions),1);
            for i=1:length(fractions)
                stats = GC.generateStats('fraction',fractions(i));
                p_thetas{i} = stats.p_theta(:);
                correlations(i) = corr(p_thetas{i},truth);
            end
            
            % convert fractions to n trials
            n_trials = round(fractions*size(GC.pf.rec,2));
            
            % convert r to r^2
            r_squared = correlations.^2;
            
            plot(n_trials,r_squared,'-o','Parent',uiParent);
            xlim([0 max(n_trials)]);
            title('Cross-Validation (r^2)');
            xlabel('trials');
            ylabel('r^2');
        end

    end

    methods(Static)
        
        function [p_theta] = findTuningStats(sts, orientations, orientation_biases)
        %FINDTUNINGSTATS compute p(theta|spike) distribution from sts
        %
        %   NB: This method is static because it is passed as a function
        %   pointer to the bootstrp() function in GridCurv.generateStats()
        %
        %   PARAM sts           (n_spike x n_grid_cell matrix) spike-triggered stimuli ensemble
        %   PARAM theta_bias    (n_grid_cell x n_ori matrix) corrects orientation anisotropy
        %   RETURN p_theta      (n_grid_cell x n_ori matrix) p(theta|spike) for each position

            % Count orientation occurrences, divide by total to get
            % p(theta|spike)
            hist_results = zeros(size(orientation_biases));
            
            N_CELLS = size(orientation_biases,1);
            for i=1:N_CELLS
                [counts,~] = hist(sts(:,i),orientations);
                hist_results(i,:) = counts;
            end
            p_theta = hist_results / sum(hist_results(1,:)); % max(sum(hist_results(1,:)),sum(hist_results(2,:)));

            % Correct for orientation bias, normalize so that sum of p at
            % each grid cell still == 1
            p_theta = p_theta .* orientation_biases;
            for i=1:N_CELLS
                 p_theta(i,:) = p_theta(i,:) * (1 / sum(p_theta(i,:)));
            end
        end
    end
    
end

