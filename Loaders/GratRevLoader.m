classdef GratRevLoader < SuperLoader
    % GRATREVLOADER loads gratrev stimuli
    
    properties
    end
    
	properties (Constant)
        PADDING_RATIO = 0.25
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
        
        function stimulus = getByStimulusParams(GRL, stimulusParams)
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
        % STIMULUS is a double in the range [0,1]
            
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
            outer_mask = (circle_xx.^2 + circle_yy.^2) < circle_radius^2;
            img(~outer_mask) = 0.5;
            
            % Add Taper
            taper_size = GRL.pf.rec(1).params.taper;
            taper_alpha = linspace(0,1,taper_size+2); % add 1 to include 0 and 1 in the sequence
            taper_alpha = taper_alpha(2:end-1);
            % @todo experiment with different taper slopes
            taper_alpha = taper_alpha.^1.5;
            
            % Alpha blend, one concentric circle at a time
            for i=1:taper_size
                circle_x = radius;
                circle_y = radius;
                circle_radius = radius - i;
                [circle_xx,circle_yy] = ndgrid((1:size(img,1))-circle_x,(1:size(img,2))-circle_y);
                % Find circle
                taper_mask = ((circle_xx.^2 + circle_yy.^2) >= circle_radius^2)...
                            & ((circle_xx.^2 + circle_yy.^2) < (circle_radius+1)^2);
                % Alpha blending
                img(taper_mask) = taper_alpha(i) * img(taper_mask) + (1-taper_alpha(i)) * 0.5;
            end
            
            % Add Padding
            padding_size = round(size(img,1)*GRL.PADDING_RATIO); % amount of padding
            bg = 0.5 * ones(size(img,1)+padding_size);
            stimulus = bg;
            pad = round(padding_size / 2);
            stimulus(pad:(pad + size(img,1) - 1),pad:(pad + size(img,1) - 1)) = img;
        end
    end
end
