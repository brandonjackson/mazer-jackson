classdef SuperLoader < handle
    % SUPERLOADER base class for all stimulus loaders
    
    properties
        % Experiment P2M File
        pf
    end
    
	properties (Constant)
        % Stimuli are loaded at this downsampled scale
        DOWNSAMPLE_SCALE = 0.5
    end
    
    methods
        
        function SL = SuperLoader(pf)
        % SUPERLOADER constructs the class
            SL.pf = pf;
        end
        
        function img = randomStimulus(SL)
        % RANDOMSTIMULUS loads a random image into memory
        %   RETURN img                  (double matrix) grayscale image
        end

    end
end