classdef GratRevLoader < SuperLoader
    % GRATREVLOADER loads gratrev stimuli
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
        function GRL = GratRevLoader(pf)
            GRL = GRL@SuperLoader(pf); % Construct superclass
        end
        
        function [img] = randomStimulus(GRL)
        % RANDOMSTIMULUS loads a random gratrev stimulus image
        % See GratRevLoader.getByStimulusParams for details of image
        % generation.
        
            % @todo load the stimulus space from the p2m file
            oris = 0:15:165;
            phases = [0,180];
            sfs = [1,2,4,8,16,32,64];
            
            stimulusParams = {};
            stimulusParams.rmult = 1;
            stimulusParams.stimulusSize = GratRevUtil.getStimulusSize(GRL.pf);
            stimulusParams.ori = oris(randi(length(oris)));
            stimulusParams.phase = phases(randi(length(phases)));
            stimulusParams.sf = sfs(randi(length(sfs)));
            stimulusParams.sf = (round(stimulusParams.sf/stimulusParams.rmult) / stimulusParams.stimulusSize);
            stimulusParams.stype = 0; % grating
            
            img = GRL.getByStimulusParams(stimulusParams);
        end
        
        function [img] = getByTrigger(GRL, pf, trigger)
            stimulusParams = GratRevUtil.trigger2stimulusParams(pf, trigger);
            img = GRL.getByStimulusParams(stimulusParams);
        end
        
        function [img] = getByStimulusParams(GRL, stimulusParams)
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
        % IMG is a double in the range [0,1]
            
            radius = stimulusParams.stimulusSize / 2; % @todo make sure this is an even number
            img = mkgrating(radius,...
                's',...
                stimulusParams.ori,...
                stimulusParams.phase,...
                stimulusParams.sf,...
                1);
            img = img / 255; % rescale from [0,255] to [0,1] range
            
            % create circular mask
            circle_x = radius;
            circle_y = radius;
            circle_radius = radius;
            [circle_xx,circle_yy] = ndgrid((1:size(img,1))-circle_x,(1:size(img,2))-circle_y);
            mask = (circle_xx.^2 + circle_yy.^2) < circle_radius^2;
            img = ((img - 0.5) .* mask + 0.5);
        end
    end
end
