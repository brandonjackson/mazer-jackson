classdef ComplexCellModel < SuperModel
    % COMPLEXCELLMODEL complex cell model
    
    properties

        ori
        
        wavelength
        
        sigma
        
        phase
        
        aspect
        
        nonlinearities
                
      %  x_offset
        
     %   y_offset
        
        SC1
        
        SC2
        
        cacheEnabled
        
        subunitInputType

    end
    
    properties (Constant)
                
    end
    
    methods

        function CX = ComplexCellModel(args, loader)
        % COMPLEXCELLMODEL constructs a new model object.
        % This constructor creates two simple cell subunits. It requires
        % multiple args, which are listed below.
        %
        % PARAMS
        %   args
        %       orientation
        %       argsgabor_params
        %       nonlinearities
        %       stimulusLoader
        %       stimulusSize
        %       x_offset
        %       y_offset
        %       imageLoader
        %       fit
        %       cacheEnabled
        %       subunitAlgorithm
       
            if nargin < 2
                loader = [];
                %warning('no stimulus loader provided');
            else
                CX.setStimulusLoader(loader);
            end
                

            % Define orientation and offsets
            CX.ori = args.ori;
            CX.wavelength = args.wavelength;
            CX.sigma = args.sigma;
            CX.phase = args.phase;
            CX.aspect = args.aspect;
            %CX.x_offset = args.x_offset;
            %CX.y_offset = args.y_offset;
            
            CX.nonlinearities = args.nonlinearities;
            CX.cacheEnabled = args.cacheEnabled;
            CX.subunitInputType = args.subunitInputType;
            
            
            
            % Simple Cell 1
            SC1_args = CX.getSimpleCellArgs('phase',0);
            CX.SC1 = SimpleCellModel(SC1_args, loader);
            
            % Save Nonlinearities For Later
            CX.nonlinearities.sc_a = CX.SC1.nonlinearities.sc_a;
            CX.nonlinearities.sc_b = CX.SC1.nonlinearities.sc_b;

            % Simple Cell 2
            SC2_args = CX.getSimpleCellArgs('phase',pi/2);
            CX.SC2 = SimpleCellModel(SC2_args, loader);
            
            if ~isfield(args,'fit') || args.fit==1
                CX.fit();
            end
        end
        
        function OUT = stimulate(CX, stimulus)

        % STIMULATE Models complex cell behavior, with caching support
        % 
        % PARAMS
        %   (optional: if none provided, uses stimulusLoader) stimulus 
        %   (optional, for SC model: response, response_sum)
        %
        % OUTPUT
        %   OUT.firing_rate
        %   OUT.spike
        %   OUT.subunit_sum
        %   OUT.SC1_firing_rate
        %   OUT.SC1_response_sum
        %   OUT.SC2_firing_rate
        %   OUT.SC2_response_sum
        
            % Load random stimulus if none provided
            if nargin < 2
                stimulus = CX.stimulusLoader.randomStimulus();
%             elseif ~isequal(size(stimulus),size(CX.stimulusSize))
%                 warning('stimulus provided does not match CX.stimulusSize');
            end

            % Output variable
            OUT = {};

            % Outline:
            % 1) Create two simple cells
            % 2) Pass them stimulus, record output firing rate
            % 3) Square Firing Rates
            % 4) Add up Firing Rates
            % 5) Rectify
            % 6) Poisson

            %% Simple Cells

            % Simple Cell One
            SC1_OUT = CX.SC1.stimulate(stimulus);
            OUT.SC1_firing_rate = SC1_OUT.firing_rate;
            OUT.SC1_response_sum = SC1_OUT.response_sum;

            % Simple Cell Two
            SC2_OUT = CX.SC2.stimulate(stimulus);
            OUT.SC2_firing_rate = SC2_OUT.firing_rate;
            OUT.SC2_response_sum = SC2_OUT.response_sum;
            
            %% Convert Subunit Responses into Membrane Potential and Firing Rates
            % The cell can respond to three different kinds of
            % `subunitInputType`s. These inputs are converted into a membrane 
            % potential (called `subunit_sum`) and then into a firing rate.
            
            % 1) `rates`: Firing Rates            
            if strcmp(CX.subunitInputType,'rates')
                % This is the "Sum" step in the "Square and Sum" algorithm
                % The squaring occurs in the simple cells, since they respond
                % to both positive and negative convolution outputs
                OUT.subunit_sum = SC1_OUT.firing_rate + SC2_OUT.firing_rate;
                
                % Convert into Firing Rate
                if ~isfield(CX.nonlinearities, 'cx_a') || ~isfield(CX.nonlinearities, 'cx_b')
                    OUT.firing_rate = OUT.subunit_sum / 2; % average subunit rates
                else
                    OUT.firing_rate = CX.nonlinearities.cx_a * OUT.subunit_sum^1.3 + CX.nonlinearities.cx_b;
                end
                
            % 2) `sums`: Filter Convolution Sums
            % This is useful when you don't want to throw away information
            % from the subunits. It reads directly from the output of the
            % linear filters.
            elseif strcmp(CX.subunitInputType,'sums')
                OUT.subunit_sum = abs(SC1_OUT.response_sum) + abs(SC2_OUT.response_sum);
                
                % Convert into Firing Rate if NL defined
                if ~isfield(CX.nonlinearities, 'cx_a') || ~isfield(CX.nonlinearities, 'cx_b')
                    OUT.firing_rate = OUT.subunit_sum;
                else
                    OUT.firing_rate = CX.nonlinearities.cx_a * OUT.subunit_sum^2 + CX.nonlinearities.cx_b;
                end
            
            % 3) `spikes`: Spikes
            else
                % @todo implement this!
            end
            
            %% Generate noisy firing rate and spike
            [OUT.spike, ~] = ComplexCellModel.rate2spike(OUT.firing_rate, CX.DT);
        end
        
        function [] = resetCache(CX)
            CX.SC1.resetCache();
            CX.SC2.resetCache();
        end
        
        function [] = fit(CX)
            
        % fit() Fit nonlinearities of a complex cell model
        %
        % PARAMS
        %   IN.IMAGE_DIR
        %   IN.image_files
        %   IN.stimulusSize
        %   IN.gabor_params
        %       kernel_size
        %       wavelength
        %       phase
        %       sigma
        %       aspect
        %       theta
        %   IN.sc_a
        %   IN.sc_b
        %
        % OUTPUTS
        %   a
        %   b
        %
        % CURRENTLY BROKEN!
            N_STIMULI = 2000;
            N_POSITIONS = 1;
            IDEAL_RATE = 40;
            
            subunit_sums = zeros(N_STIMULI*N_POSITIONS,1);
            
            i = 1;

            for s=1:N_STIMULI

                stim = CX.stimulusLoader.randomStimulus();
                stim = stim - mean(stim(:));

                for pos=1:N_POSITIONS
                    xmin = (round(size(stim,2)/2.7)) - (CX.gabor_params.kernel_size/2);
                    
                    CX_args = CX.getSimpleCellArgs();
                    CX_args.fit = 0;
                    CX_args.orientation = rand()*pi;
                    CX_args.subunitInputType = CX.subunitInputType;
                    
                    if pos==1
                        CX_args.x_offset = 0;
                        CX_args.y_offset = 0;
                    else
                        CX_args.x_offset = randi([-1*xmin xmin]);
                        CX_args.y_offset = randi([-1*xmin xmin]);
                    end
                    
                    CX_temp = ComplexCellModel(CX_args);
                    CX_OUT = CX_temp.stimulate(stim);
                    subunit_sums(i) = CX_OUT.subunit_sum;
                    i = i + 1;
                end
            end
            min_threshold = max(abs(subunit_sums))*0.03;
            fprintf('[cx_fit] n subunit_sums (pre-screening) = %f\n',length(subunit_sums));
            subunit_sums = subunit_sums(abs(subunit_sums) > min_threshold);
            fprintf('[cx_fit] n subunit_sums (post-screening) = %f\n',length(subunit_sums));

            ideal = fminsearch(@prod2rate,[1.0 0.0],[], subunit_sums, 1.3, IDEAL_RATE);
            CX.nonlinearities.cx_a = ideal(1);
            CX.nonlinearities.cx_b = ideal(2);

            figure();
            subplot(121);
            hist(subunit_sums(:),40);
            title('CX Subunit FR Sums');
            subplot(122);
            fr = CX.nonlinearities.cx_a*subunit_sums.^1.3 + CX.nonlinearities.cx_b;
            fprintf('[cx_fit] FR Stats (pre-ceiling): mean=%.0f, median=%.0f\n',mean(fr(:)),median(fr(:)));
            fr(fr > 200) = 200;
            fprintf('[cx_fit] FR Stats (post-ceiling): mean=%.0f, median=%.0f\n',mean(fr(:)),median(fr(:)));
            hist(fr(:),40);
            title('CX Firing Rates');
        end
        
        function [SC_args] = getSimpleCellArgs(CX, varargin)
            SC_args = {};
            SC_args.ori = CX.ori;
            SC_args.wavelength = CX.wavelength;
            SC_args.sigma = CX.sigma;
            SC_args.phase = CX.phase;
            SC_args.aspect = CX.aspect;
            %SC_args.x_offset = CX.x_offset;
            %SC_args.y_offset = CX.y_offset;
            SC_args.stimulusLoader = CX.stimulusLoader;
            SC_args.stimulusSize = CX.stimulusSize;
            SC_args.nonlinearities = CX.nonlinearities;   
            SC_args.fullwaverectify = 1;
            SC_args.cacheEnabled = CX.cacheEnabled;
            
            if ~isempty(varargin) && strcmp(varargin{1},'phase') && length(varargin) > 1
                SC_args.phase = varargin{2};
            end
        end

        function [] = plotKernel(CX)
            figure();
            subplot(1,2,1);
            imagesc(CX.SC1.kernel);
            axis square;
            title('SC1 Kernel');
            colorbar();
            caxis([-1 1]);
            subplot(1,2,2);
            imagesc(CX.SC2.kernel);
            axis square;
            colorbar();
            caxis([-1 1]);
            title('SC2 Kernel');
            boxtitle(sprintf('CX Kernel \\lambda = %dpx, \\sigma = %dpx', CX.wavelength, CX.sigma));
        end
    end
    
    methods (Static)
        
        function [] = plotQuickKernel(ori, wavelength, sigma)
            
            gb1 = GaborUtil.createFilter(ori, wavelength, sigma, 0);
            gb2 = GaborUtil.createFilter(ori, wavelength, sigma, pi/2);
            
            figure();
            subplot(1,2,1);
            imagesc(gb1);
            axis square;
            title('SC1 Kernel');
            colorbar();
            caxis([-1 1]);
            subplot(1,2,2);
            imagesc(gb2);
            axis square;
            colorbar();
            caxis([-1 1]);
            title('SC2 Kernel');
            boxtitle(sprintf('CX Kernel \\lambda = %dpx, \\sigma = %dpx', wavelength, sigma));
            
        end
        
    end
end