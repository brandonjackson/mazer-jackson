classdef GratRevBarLoader < SuperLoader
    % GRATREVBARLOADER loads gratrev bar stimuli
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
        function GRBL = GratRevBarLoader(pf)
           GRBL = GRBL@SuperLoader(pf); % Construct superclass
        end
        
        function [img] = randomStimulus(GRBL)
        % RANDOMSTIMULUS loads a random gratrev stimulus image
        % See GratRevLoader.getByStimulusParams for details of image
        % generation.
        
            % @todo load the stimulus space from the p2m file
            oris = 0:15:165;
            phases = [0,180];
            sfs = [1,2,4,8,16,32,64]; % nb: width in pixels
            
            stimulusParams = {};
            stimulusParams.rmult = 1;
            stimulusParams.stimulusSize = GratRevUtil.getStimulusSize(GRL.pf);
            stimulusParams.ori = oris(randi(length(oris)));
            stimulusParams.phase = phases(randi(length(phases)));
            stimulusParams.sf = sfs(randi(length(sfs)));
            stimulusParams.stype = 1; % bar
            
            img = GRBL.getByStimulusParams(stimulusParams);
        end
        
        function [img] = getByStimulusParams(GRBL, stimulusParams)
        % GETBYSTIMULUSPARAMS loads a bar image based on the details in
        % STIMULUSPARAMS. The STIMULUSPARAMS struct can be derived from the
        % ev_e triggers using the GratRevUtil.trigger2stimulusParams()
        % function. 
        %
        % STIMULUSPARAMS contains the following properties:
        %   - stimulusSize      (int) width of square stimulus
        %   - ori               (double) orientation
        %   - phase             (double) phase
        %   - sf                (double) bar width, in px
        %
        % IMG is a double in the range [0,1]
            img = zeros(stimulusParams.stimulusSize,'double');
            width = stimulusParams.sf;
            center = stimulusParams.stimulusSize / 2;
            if stimulusParams.phase==0
                fill = -0.5;
            else
                fill = 0.5;
            end
            
            % Draw rectangle in left corner
            img(:,1:width) = fill;
            
            % Then shift to center
            img = circshift(img,[0 round(center - (width/2))]);
            
            % Rotate
            img = imrotate(img,stimulusParams.ori,'bilinear','crop');
            
            % Add 0.5, so that image range is [0,1]
            % Note: the range was [-0.5,0.5] before so that when the image
            % was rotated the pixels added in the corner had the background
            % color of 0
            img = img + 0.5; 
            
            % @todo add gaussian envelope?
        end
    end
end