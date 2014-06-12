classdef SimpleCellModel < SuperModel
    % SIMPLECELLMODEL complex cell model
    
    properties

        orientation
        
        phase
        
        gabor_params
        
        fullwaverectify
        
        nonlinearities
        
        kernel
        
        kernel_pos
        
        x_offset
        
        y_offset
        
        cacheEnabled
                
        cached_response
        
        cached_response_sum
        
        cached_firing_rate

    end
    
    properties (Constant)
        
    end
    
    methods
        function [ SC ] = SimpleCellModel(args)
        % SIMPLECELLMODEL constructs a new model object.
        % This constructor creates two simple cell subunits. It requires
        % multiple args, which are listed below.
        %
        % PARAMS
        %   args
        %       orientation
        %       phase
        %       gabor_params
        %       nonlinearities
        %       stimulusLoader
        %       stimulusSize
        %       x_offset
        %       y_offset

            SC.orientation = args.orientation;
            SC.phase = args.phase;
            SC.gabor_params = args.gabor_params;
            SC.nonlinearities = args.nonlinearities;
            SC.setStimulusLoader(args.stimulusLoader);
            SC.x_offset = args.x_offset;
            SC.y_offset = args.y_offset;

            if ~isfield(args,'fullwaverectify')
                SC.fullwaverectify = 0;
            else
                SC.fullwaverectify = args.fullwaverectify;
            end
            
            if ~isfield(args,'cacheEnabled')
                SC.cacheEnabled = 0;
            else
                SC.cacheEnabled = args.cacheEnabled;
            end

            % Create Filter
            gabor = gabor_filter(args.gabor_params,'phase',args.phase,'theta',args.orientation);
            
            [SC.kernel,SC.kernel_pos] = gabor_filter_translate(gabor, args.stimulusSize, args.x_offset, args.y_offset);
            
            if SC.nonlinearities.sc_a==1 && SC.nonlinearities.sc_b==0
                SC.fit();
            end
        end
        
        function [ OUT ] = stimulate( SC, stimulus )
        % STIMULATE Models simple cell behavior, with caching support
        % 
        % PARAMS
        %   stimulus
        %
        % OUTPUT
        %   OUT.firing_rate
        %   OUT.spike
        %   OUT.response
        %   OUT.response_sum
        %   OUT.final_kernel

            % Poisson Scale
            POISSON_SCALE = 0.001;

            % Output variable
            OUT = {};
          %  SC.stimulusSize
            SC.stimulusSize = size(stimulus);

            % Grab stimulus frame and convolve with kernel
            if(~SC.hasCache())
                kernel_size = size(SC.kernel);

                if isequal(SC.stimulusSize,kernel_size)
                    kernel = SC.kernel;
                
                % stimulus smaller than kernel
                elseif SC.stimulusSize(1) < kernel_size(1)
                    pad = round((kernel_size(1) - SC.stimulusSize(1)) / 2);
                    stimulus_padded = zeros(kernel_size) + 0.5;
                    stimulus_padded((pad+1):(SC.stimulusSize(1)+pad),(pad+1):(SC.stimulusSize(1)+pad)) = stimulus;
                    stimulus = stimulus_padded;
                    kernel = SC.kernel;
                
                % stimulus bigger than kernel
                else
                    origin = SC.stimulusSize/2;
                    left = origin(2) - floor(kernel_size(2)/2);
                    right = origin(2) + floor(kernel_size(2)/2);
                    top = origin(1) - floor(kernel_size(1)/2);
                    bottom = origin(1) + floor(kernel_size(1)/2);
                    kernel = zeros(SC.stimulusSize);
                    kernel(left:right,top:bottom) = SC.kernel;
                end

                OUT.final_kernel = kernel;

                % only convolve with the relevant slice of the stimulus
                % kernel_pos is [x1 x2 y1 y2] of slice
                x_slice = SC.kernel_pos(1):SC.kernel_pos(2);
                y_slice = SC.kernel_pos(3):SC.kernel_pos(4);
                kernel_slice = kernel(x_slice,y_slice);
                stimulus_slice = stimulus(x_slice,y_slice);
                OUT.response = kernel_slice .* double(stimulus_slice);

                % Save Sum
                OUT.response_sum = sum(OUT.response(:));

                % Rectify, Generate Firing Rates
                % Full-wave rectification makes resposne phase invariant
                if(SC.fullwaverectify==1)
                    rectified_response_sum = abs(OUT.response_sum);
                else
                    rectified_response_sum = OUT.response_sum;
                end

                OUT.firing_rate = SC.nonlinearities.sc_a * OUT.response_sum^2 + SC.nonlinearities.sc_b;
                
                if SC.cacheEnabled==1
                    SC.cached_firing_rate = OUT.firing_rate;
                    SC.cached_response_sum = OUT.response_sum;
                    SC.cached_response = OUT.response;
                end

            % Use Cached Responses
            else
                OUT.response = SC.cached_response;
                OUT.response_sum = SC.cached_response_sum;
                OUT.firing_rate = SC.cached_firing_rate;
            end

            % Add Noise, Generate Spike From Firing Rate
            [OUT.spike, OUT.firing_rate] = SC.rate2spike(OUT.firing_rate, SC.DT);
        end
        
        function [result] = hasCache(SC)
            result = ~isempty(SC.cached_firing_rate) && ~isempty(SC.cached_response) && ~isempty(SC.cached_response_sum);
        end
        
        function [] = resetCache(SC)
            SC.cached_firing_rate = [];
            SC.cached_response = [];
            SC.cached_response_sum = [];
        end
        
        function [] = fit(SC)
        % fit() Fits nonlinearity for a simple cell
        % saves to:
        %   SC.nonlinearities.sc_a
        %   SC.nonlinearities.sc_b
        %
        % PARAMS
        %   image_dir

            N_STIMULI = 2000;
            N_POSITIONS = 1;
            IDEAL_RATE = 20;
            
            % Create Filter
            

            response_sums = zeros(N_STIMULI * N_POSITIONS,1);

            i = 1;

            for s=1:N_STIMULI

                stim = SC.stimulusLoader.randomStimulus();
                stim = stim - mean(stim(:));
                SC.gabor_params.theta = rand() * pi;
                gb = gabor_filter(SC.gabor_params);

                for pos=1:N_POSITIONS
                    xmin = (round(size(stim,2)/2.7)) - (SC.gabor_params.kernel_size/2);
                    
                    % start at center
                    if pos==1
                        x_offset = 0;
                        y_offset = 0;
                    else
                        x_offset = randi([-1*xmin xmin]);
                        y_offset = randi([-1*xmin xmin]);
                    end
                    
                    new_gb = gabor_filter_translate(gb, size(stim), x_offset ,y_offset);

                    response = double(stim) .* new_gb;

                    response_sums(i) = sum(response(:));
                    i = i + 1;
                end
            end
            
            min_threshold = max(abs(response_sums))*0.03;
            fprintf('[sc_fit] n response_sums (pre-screening) = %f\n',length(response_sums));
            response_sums = response_sums(abs(response_sums) > min_threshold);
            fprintf('[sc_fit] n response_sums (post-screening) = %f\n',length(response_sums));

            ideal = fminsearch(@prod2rate,[1.0 0.0],[], response_sums, 2, IDEAL_RATE);
            SC.nonlinearities.sc_a = ideal(1);
            SC.nonlinearities.sc_b = ideal(2);

            figure();
            subplot(121);
            hist(response_sums(:),40);
            title('SC Response Sums');
            subplot(122);
            fr = SC.nonlinearities.sc_a*response_sums.^2 + SC.nonlinearities.sc_b;
            %fr(fr > 200) = 200;
            hist(fr(:),20);
            
            title(['SC Firing Rates (mean=' num2str(round(mean(fr))) ')']);
        end

    end
    
    methods(Static)
        
    end
end