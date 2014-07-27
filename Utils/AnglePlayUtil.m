classdef AnglePlayUtil < SuperUtil
    % ANGLEPLAYUTIL library of static methods for handling curvplay p2m files
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
    end
    
    methods (Static)
        
%         function [stimuli, responses] = getResponses(pf, latency, winsize)
%         % GETRESPONSES extracts a STIMULI and RESPONSES vector for a given
%         % p2m file PF using the LATENCY and WINSIZE. The STIMULI vector
%         % contains a stimulusParams struct that describes each stimulus,
%         % and the RESPONSES vector contains the firing rate in spikes/s.
%         
%             rast = prast(pf,'pre',latency,'post',latency + winsize);
%             
%             % Stimuli vector
%             stimuli = cell(length(rast.triggers),1);
%             for i=1:length(rast.triggers)
%                 stimuli{i} = GratRevUtil.trigger2stimulusParams(pf, rast.triggers{i});
%             end
%             
%             % Responses Vector (in spikes/s)
%             spike_counts = nansum(rast.data,2);
%             responses = spike_counts / (winsize * GratRevUtil.DT);
%         end
        
%         function stimulusSize = getStimulusSize(pf)
%             stimulusSize = (2 * pf.rec(1).params.radius + pf.rec(1).params.taper);
%         end
        
        function stimulusLoader = getStimulusLoader(pf)
        % GETSTIMULUSLOADER set up stimulus loader of type bar or grating
        % Works by looking at the stimulus type parameter in the pf file
            stimulusLoader = AnglePlayLoader();
        end
        
        function version = taskVersion(pf)
        % TASKVERSION determines whether old or new version of task was
        % used
        %
        % The function grabs the revision of the tasks repo that was used
        % when the pf was recorded.
        %
        % Output:
        %   version - 0 if old version, 1 if new version
        
            version_str = strjoin(pf.rec(1).params.X_version_info);
            version_match = regexp(version_str,'curvplay.py\s(\d*)','tokens');
            revision_number = str2double(strjoin(version_match{:}));
            version = ~(revision_number < 938);
        end
        
        function stimulusParams =  trigger2stimulusParams(pf, trigger)
        %% Extract the stimuli associated with each ev_e trigger
        % Input: 
        % - pf: struct with p2m file data
        % - trigger: string of ev_e trigger
        % Output:
        % - stimulusParams: struct with 6 stimulus parameters:
        % { ori, sf, phase, contrast, rmult (radius multiplier), stype (stimulus type)}

            inum = sscanf(trigger, '%*s %d', 1);

            % all info in IMAGE_INFO indexed by inum, type(abs or
            % parabola),coefficient,stroke,rotation,ori,scale,polarity,image_file_n
            % e.g. 'image-000006.png parabola 2.000000 7.000000 -90.0 666.0 1.00 1 6'
            info = pf.rec(1).params.IMAGE_INFO{inum+1}; % switch from 0- to 1-based
            
            tokens = textscan(info,'%s %s %f %f %f %f %f %d %d');
            %tokens = strsplit(info,[char(9) ' ']); % split by tab+space

            % Extract and accumulate stimulus parameters            
            stimulusParams = {};

            % Parabola == 0, Abs == 1
            type_label = tokens{2};
            if strcmp(type_label,'parabola')
                stimulusParams.type = 0;
            else
                stimulusParams.type = 1;
            end

            % note: ignore the EV_ORI param and compute it from
            % the rotation value -- base image is horiz (0deg)
            % and negative rotation values indicate CCW rotation
            ev_rotation = tokens{5};
            ev_ori = tokens{6};
            stimulusParams.ori = 90 - ev_rotation - 360;

            stimulusParams.coefficient = tokens{3};
            stimulusParams.stroke = tokens{4};
            stimulusParams.scale = tokens{7};
            stimulusParams.polarity = tokens{8};
            stimulusParams.image_file_n = tokens{9};
            stimulusParams.inum = inum;
            %stimulusParams.stimulusSize = AnglePlayUtil.getStimulusSize(pf);
            
            % Compute Angle
            if stimulusParams.type==1 % abs
                 stimulusParams.angle = 2*atand(1/stimulusParams.coefficient);
            else % parabola: find where they intercept a circle with r=2
                x = max(real(roots([stimulusParams.coefficient^2 0 1 0 -4])));
                y = stimulusParams.coefficient*(x^2);
                stimulusParams.angle = 2*atand(x/y);
            end
        end
    end
end