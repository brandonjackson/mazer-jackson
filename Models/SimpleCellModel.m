classdef SimpleCellModel < SuperModel
    % SIMPLECELLMODEL complex cell model
    
    properties

        ori
        
        wavelength
        
        sigma
        
        phase
        
        aspect
                
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
        function [ SC ] = SimpleCellModel(args, loader)
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
        
            if nargin < 2 || isempty(loader)
                loader = [];
               %warning('no stimulus loader provided');
            else
                SC.setStimulusLoader(loader);
            end

            SC.ori = args.ori;
            SC.wavelength = args.wavelength;
            SC.sigma = args.sigma;
            SC.phase = args.phase;
            SC.aspect = args.aspect;
            SC.nonlinearities = args.nonlinearities;
            SC.x_offset = 0;
            SC.y_offset = 0;
            
            

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
            gabor = GaborUtil.createFilter(args.ori, args.wavelength, args.sigma, args.phase, args.aspect);
            %gabor = gabor_filter(args.gabor_params,'phase',args.phase,'theta',args.orientation);
            SC.kernel = gabor;
            
            %[SC.kernel,SC.kernel_pos] = gabor_filter_translate(gabor, args.stimulusSize, args.x_offset, args.y_offset);
            
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
            
            SC.stimulusSize = size(stimulus);

            % Grab stimulus frame and convolve with kernel
            if(~SC.hasCache())
                kernel_size = size(SC.kernel);

                if isequal(SC.stimulusSize,kernel_size)
                    OUT.final_kernel = SC.kernel;
                
                % Stimulus smaller than kernel: cut out slice of kernel
                % that is the size of the stimulus (since everything
                % outside of the stimulus would be zero if we added padding
                % to the stimulus anyways)
                elseif SC.stimulusSize(1) < kernel_size(1)
                    margin = round((kernel_size - SC.stimulusSize(1)) / 2);
                    OUT.final_kernel = SC.kernel(margin:(SC.stimulusSize(1)+margin-1),margin:(SC.stimulusSize(2)+margin-1));
                
                % Stimulus bigger than kernel: cut out slice of stimulus
                % that is the size of the kernel
                else
                    margin = round((SC.stimulusSize(1) - kernel_size) / 2);
                    stimulus = stimulus(margin:(kernel_size(1)+margin-1),margin:(kernel_size(2)+margin-1));
                    OUT.final_kernel = SC.kernel;
                end
                
                stimulus = stimulus - 0.5; % center

                OUT.response = OUT.final_kernel .* stimulus;

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
        %
        % CURRENTLY BROKEN!

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