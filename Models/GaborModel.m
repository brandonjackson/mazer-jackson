classdef GaborModel < SuperModel
    % GABORMODEL model of synthetic orientation-selective cell with gabor kernel
    
    properties

        orientation
        
        gabor_params
        
        nonlinearities
        
        x_offset
        
        y_offset
        
    end
    
    methods

        function GM = GaborModel(args)
        % GABORMODEL constructs a new model object.
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
        %       fit

            % Define orientation and offsets
            GB.orientation = args.orientation;
            GB.gabor_params = args.gabor_params;
            GB.x_offset = args.x_offset;
            GB.y_offset = args.y_offset;
            GB.setStimulusLoader(args.stimulusLoader);
            GB.nonlinearities = args.nonlinearities;

            % Simple Cell 1
            SC1_args = GB.getSimpleCellArgs('phase',0);
            GB.SC1 = SimpleCellModel(SC1_args);
            
            % Save Nonlinearities For Later
            GB.nonlinearities.sc_a = GB.SC1.nonlinearities.sc_a;
            GB.nonlinearities.sc_b = GB.SC1.nonlinearities.sc_b;

            % Simple Cell 2
            SC2_args = GB.getSimpleCellArgs('phase',pi/2);
            GB.SC2 = SimpleCellModel(SC2_args);
            
            if ~isfield(args,'fit') || args.fit==1
                GB.fit();
            end
        end
        
        function OUT = stimulate(CX, stimulus)

        % STIMULATE Models complex cell behavior, with caching support
        % 
        % PARAMS
        %   stimulus
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

            % Poisson Scale (usually == dt)
            POISSON_SCALE = 0.001;

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

            %% Square and Sum
            OUT.subunit_sum = SC1_OUT.firing_rate + SC2_OUT.firing_rate;

            if ~isfield(CX.nonlinearities, 'cx_a') || ~isfield(CX.nonlinearities, 'cx_b')
                OUT.firing_rate = OUT.subunit_sum / 2;
            else
                OUT.firing_rate = CX.nonlinearities.cx_a * OUT.subunit_sum^1.3 + CX.nonlinearities.cx_b;
            end

            % Add Noise to Firing Rate
            OUT.firing_rate = OUT.firing_rate * (1 + 0.1*randn());

            % Generate Spike via Poisson Process
            if rand() < POISSON_SCALE * OUT.firing_rate
                OUT.spike = 1;
            else
                OUT.spike = 0;
            end
        end
        
        function [] = fit(CX)
            
        % cx_fit() Fit nonlinearities of a complex cell model
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

    end
    
    methods(Static)
        
        function filtered_img = filter2_padded(h, img)
            img_width = size(img,1);
            padding = round(img_width * 0.33);

            img_padded = ones(img_width + (2*padding)) * max(img(:));
            img_padded(padding:(img_width+padding-1), padding:(img_width+padding-1)) = img;

            filtered_img_padded = filter2(h, img_padded);
            filtered_img = filtered_img_padded(padding:(img_width+padding-1), padding:(img_width+padding-1));
        end
        
    end
end