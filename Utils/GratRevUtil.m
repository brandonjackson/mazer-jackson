classdef GratRevUtil < SuperUtil
    % GRATREVUTIL library of static methods for handling gratrev p2m files
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
    end
    
    methods (Static)
        
        function [stimuli, responses] = getResponses(pf, latency, winsize)
        % GETRESPONSES extracts a STIMULI and RESPONSES vector for a given
        % p2m file PF using the LATENCY and WINSIZE. The STIMULI vector
        % contains a stimulusParams struct that describes each stimulus,
        % and the RESPONSES vector contains the firing rate in spikes/s.
        
            rast = prast(pf,'pre',latency,'post',latency + winsize);
            
            % Stimuli vector
            stimuli = cell(length(rast.triggers),1);
            for i=1:length(rast.triggers)
                stimuli{i} = GratRevUtil.trigger2stimulusParams(pf, rast.triggers{i});
            end
            
            % Responses Vector (in spikes/s)
            spike_counts = nansum(rast.data,2);
            responses = spike_counts / (winsize * GratRevUtil.DT);
        end
        
        function stimulusSize = getStimulusSize(pf)
            stimulusSize = (2 * pf.rec(1).params.radius + pf.rec(1).params.taper);
        end
        
        function stimulusLoader = getStimulusLoader(pf)
        % GETSTIMULUSLOADER set up stimulus loader of type bar or grating
        % Works by looking at the stimulus type parameter in the pf file
            if GratRevUtil.isGrating(pf)
                stimulusLoader = GratRevLoader();
            else
                stimulusLoader = GratRevBarLoader();
            end
        end
        
        function stimulusParams =  trigger2stimulusParams(pf, trigger)
        %% Extract the stimuli associated with each ev_e trigger
        % Input: 
        % - pf: struct with p2m file data
        % - trigger: string of ev_e trigger
        % Output:
        % - stimulusParams: struct with 6 stimulus parameters:
        % { ori, sf, phase, contrast, rmult (radius multiplier), stype (stimulus type)}


            %% Check to see if this is really a gratrev/ncgratrev file:
            if ~isempty(strfind(pf.rec(1).params.X_version_info{2}, 'ncgratrev'))
              nc = 1;
            elseif ~isempty(strfind(pf.rec(1).params.X_version_info{2}, 'gratrev'))
              nc = 0;
            else
              error('pgratrev:NotGratrev', 'Not a gratrev or ncgratrev file: %s\n', ...
                    pf.rec(1).params.X_version_info{2});
            end

            %% Extract and accumulate stimulus parameters            
            stimulusParams = {};
            stimulusParams.stimulusSize = GratRevUtil.getStimulusSize(pf);

            if(nc)
                % NCGratrev has both cartesian and non-cartesian stimuli. The cartesian
                % stimuli are of the form 'FLIP s <ori> <sf>' (the task code fixes phase,
                % contrast, and rmult to 0, 1, and 1, respectively). 
                stimulusParams.phase = 0.0; 
                stimulusParams.contrast = 1.0; 
                stimulusParams.rmult = 1.0; 
                stimulusParams.stype = 0;

                e = nthword(trigger);
                if e{2}  == 's'
                    e = str2num(char(e{3:end})); %#ok<ST2NM>
                    stimulusParams.ori = e(1); 
                    stimulusParams.sf = e(2) / stimulusParams.stimulusSize;
                else
                    % Noncartesian stim are totally ignored. Use ncgratrev instead.
                    % continue; % @todo make sure this isn't a problem
                    % since we can't break since we're not in a for loop
                    % like when jamie wrote it originally
                end
            else
                % Gratrev encodes are 'FLIP <ori> <sf> <phase> <contrast>
                % <rmult>'. Not all gratrev files have contrast (from May 2007),
                % rmult (from March 2011), or stype (May 2011) entries.

                % Additionally, the gratrev 'type' (0=sine,1=bar,2=sq. wave
                % grating) affects how the spatial frequency term is interpreted.
                if ~isfield(pf.rec(1).params, 'type')
                    stimulusParams.stype = 0;
                elseif ~isfield(pf.rec(1).params, 'X_type')
                    stimulusParams.stype = pf.rec(1).params.type;
                else
                    stimulusParams.stype = -1; %This *must* be updated below.
                end

                stimulusParams.contrast = 1.0;
                stimulusParams.rmult = 1.0;

                e = nthword(trigger);
                e = str2num(char(e{2:end})); %#ok<ST2NM>
                
                stimulusParams.ori = e(1); 
                stimulusParams.sf = e(2); 
                stimulusParams.phase = e(3);

                if length(e) >= 4 
                    stimulusParams.contrast = e(4);
                end

                if length(e) >= 5 
                    stimulusParams.rmult = e(5);        
                end

                if length(e) >=6
                    stimulusParams.stype = e(6);
                end

                if stimulusParams.stype == 1
                    %Size of bar (pixels?) 
                    stimulusParams.sf = round(stimulusParams.stimulusSize/round(stimulusParams.sf/stimulusParams.rmult));
                else
                    % cycles/stim -> cycles/pixel for sine/sq. wave gratings. 
                    % NOTE: The "right" way to do this is to recalcuate stimulusParams.stimulusSize
                    % (2*radius*rmult) + taper. HOWEVER, this gets you into
                    % floating point hell, so we just divide sf by rmult, which is
                    % actually reverses what the python task does. This seems to
                    % avoid fp drama since the first division should yield an int
                    stimulusParams.sf = (round(stimulusParams.sf/stimulusParams.rmult) / stimulusParams.stimulusSize);
                end
            end
        end
        
        function stype = getStimulusType(pf)
            if ~isfield(pf.rec(1).params, 'type')
                stype = 0;
            elseif ~isfield(pf.rec(1).params, 'X_type')
                stype = pf.rec(1).params.type;
            else
                stype = pf.rec(1).params.X_type;
                %error('could not find type');
            end
        end
        
        function truth = isBar(pf)
            truth = GratRevUtil.getStimulusType(pf)==1;
        end
        
        function truth = isGrating(pf)
            truth = GratRevUtil.getStimulusType(pf)==0;
        end        
        
    end
end