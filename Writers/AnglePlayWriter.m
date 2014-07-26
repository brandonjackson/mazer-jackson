classdef AnglePlayWriter < SuperWriter
    % ANGLEPLAYWRITER writes spike times to p2m files. This modifies the old
    % file, often using predicted spikes from custom models.    
    
    properties
        % Inherited from SuperWriter:
        % - stimulusLoader
        % - pf
    end
    
	properties (Constant)
        ANGLEPLAY_PATH = '/lab/stimuli/curvplay_jackson/angles_medium-0.25';
        LATENCY = 70;
    end
    
    methods
        function APW = AnglePlayWriter(pf)
        % ANGLEPLAYWRITER constructs an object that writes to PF, a p2m file.
        %
        % PARAMS
        %    pf

            APW.pf = pf;
            APW.stimulusLoader = AnglePlayLoader();
        end
        
        function new_pf = writeFromModel(APW, model)
        % WRITEFROMMODEL generates a new p2m file, NEW_PF, using the stimuli from
        % APW.pf and a model
        %
        % NB: The model must have a method called stimulate(), which returns a
        % struct with a property called spike which is 0 or 1.
            
            new_pf = APW.pf;
            
            % Outline:
            % 1. Iterate over each PF.rec
            % 2. Erase its spike_times
            % 3. Iterate over each frameflip event
            % 4. Load stimulus image
            % 5. Iterate ms by ms until next gap_frameflip evenet. If 
            %    spike occurs, add current time + LATENCY to spike_times
            
            for i=1:length(new_pf.rec)
                new_pf.rec(i).spike_times = [];
                
                [flip_ev_i, flip_t] = p2mFindEvents(new_pf,i,'frameflip');
                [~, gap_t] = p2mFindEvents(new_pf,i,'gap_frameflip');
                
                % Iterate over each flip
                for flip_i=1:length(flip_ev_i)
                    
                    if length(gap_t) < flip_i
                        continue;
                    end
                    
                    % Calculate time interval between flip and gap flip
                    duration = gap_t(flip_i) - flip_t(flip_i);
                    
                    % Load Stimulus                    
                    ev_e = new_pf.rec(i).ev_e{flip_ev_i(flip_i)};
                    stimulus_inum = sscanf(ev_e,'frameflip %d');
                    stimulus = APW.stimulusLoader.getByImageNumber(stimulus_inum);
                    
                    % Iterate ms by ms, generating spikes
                    for t=1:duration
                        response = model.stimulate(stimulus);
                        if response.spike==1
                            spike_time = flip_t(flip_i) + (t - 1) + APW.LATENCY;
                            new_pf.rec(i).spike_times = [new_pf.rec(i).spike_times; spike_time];
                        end
                    end
                end
            end
        end
    end
end