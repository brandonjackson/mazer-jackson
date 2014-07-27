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
        
        function [] = plotRandomStimuli(SL)
        % PLOTRANDOMSTIMULI displays a grid of 12 random stimuli
            figure();
            for i=1:12
                subplot(3,4,i);
                stim = SL.randomStimulus();
                imshow(stim);
                axis square;
            end
            boxtitle(sprintf('%s Random Stimuli (%d x %d px)',...
                PFUtil.experName(SL.pf),...
                size(stim,1),...
                size(stim,2)));
        end

    end
end