classdef GratRevLoader < handle
    % GRATREVLOADER loads gratrev stimuli
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
        function GRL = GratRevLoader()
           
        end
        
        function [img] = randomStimulus(GRL, downsampled_size)
        % RANDOMSTIMULUS loads a random gratrev stimulus image
        % See GratRevLoader.getByStimulusParams for details of image
        % generation.
        
            oris = 0:15:165;
            phases = [0,180];
            sfs = [1,2,4,8,16,32,64];
            
            stimulusParams = {};
            stimulusParams.rmult = 1;
            stimulusParams.stimulusSize = 164;
            stimulusParams.ori = oris(randi(length(oris)));
            stimulusParams.phase = phases(randi(length(phases)));
            stimulusParams.sf = sfs(randi(length(sfs)));
            stimulusParams.sf = (round(stimulusParams.sf/stimulusParams.rmult) / stimulusParams.stimulusSize);
            stimulusParams.stype = 0; % grating
            
            if nargin < 2
                img = GRL.getByStimulusParams(stimulusParams);
            else
                img = GRL.getByStimulusParams(stimulusParams, downsampled_size);
            end
        end
        
        function [img] = getByTrigger(GRL, pf, trigger, downsampled_size)
            stimulusParams = GratRevUtil.trigger2stimulusParams(pf, trigger);
            if nargin < 4
                img = GRL.getByStimulusParams(stimulusParams);
            else
                img = GRL.getByStimulusParams(stimulusParams, downsampled_size);
            end
        end
        
        function [img] = getByStimulusParams(GRL, stimulusParams, downsampled_size)
        % GETBYSTIMULUSPARAMS loads a grating image based on the details in
        % STIMULUSPARAMS. The STIMULUSPARAMS struct can be derived from the
        % ev_e triggers using the GratRevUtil.trigger2stimulusParams()
        % function. 
        %
        % STIMULUSPARAMS contains the following properties:
        %   - stimulusSize      (int) width of square stimulus
        %   - ori               (double) orientation
        %   - phase             (double) phase
        %   - sf                (double) spatial frequency, in cycles/pixel
        %
        % IMG is a double in the range [-0.5, 0.5]
            
            radius = stimulusParams.stimulusSize / 2; % @todo make sure this is an even number
            img = mkgrating(radius,...
                's',...
                stimulusParams.ori,...
                stimulusParams.phase,...
                stimulusParams.sf,...
                1);
            img = im2double(img); % rescale from [0,255] to [0,1] range
            img = img - 0.5; % effectively demean
            
            if nargin < 3
                % create circular mask
                circle_x = radius;
                circle_y = radius;
                circle_radius = radius;
                [circle_xx,circle_yy] = ndgrid((1:size(img,1))-circle_x,(1:size(img,2))-circle_y);
                mask = (circle_xx.^2 + circle_yy.^2) < circle_radius^2;
                img = img .* mask;
            
            % Downsample image
            else
                img = imresize(img,[downsampled_size, downsampled_size]);
            end
        end
    end
end