classdef GratRevWriter < SuperWriter
    % GRATREVWRITER writes spike times to p2m files. This modifies the old
    % file, often using predicted spikes from custom models.
    
    properties
        % Inherited from SuperWriter:
        % - stimulusLoader
        % - pf
        % - stimulusSize
    end
    
	properties (Constant)
        LATENCY = 70;
    end
    
    methods
        function GRW = GratRevWriter(pf)
        % GRATREVWRITER constructs an object that writes to PF, a p2m file.
        %
        % PARAMS
        %    pf

            GRW.pf = pf;
            GRW.stimulusLoader = GratRevUtil.getStimulusLoader(GRW.pf);
            GRW.stimulusSize = GratRevUtil.getStimulusSize(GRW.pf);
        end
        
        function new_pf = writeFromModel(GRW, model)
        % WRITEFROMMODEL generates a new p2m file, NEW_PF, using the stimuli from
        % GRW.pf and a MODEL
        %
        % NB: The model must have a method called stimulate(), which returns a
        % struct with a property called spike which is 0 or 1.
          
            new_pf = GRW.pf;
            
            % Outline:
            % 1. Iterate over each PF.rec
            % 2. Erase its spike_times
            % 3. Iterate over each frameflip event
            % 4. Load stimulus image
            % 5. Iterate ms by ms until next gap_frameflip evenet. If 
            %    spike occurs, add current time + LATENCY to spike_times
            
            for i=1:length(new_pf.rec)
                new_pf.rec(i).spike_times = [];
                
                [flip_ev_i, flip_t] = p2mFindEvents(new_pf,i,'FLIP');
                
                % Iterate over each flip
                for flip_i=1:length(flip_ev_i)
                    
                    % Calculate time interval between flips
                    if flip_i == length(flip_ev_i)
                        % if its the last flip, use default duration
                        % @todo calculate intelligently
                        duration = 100; 
                    else
                        duration = flip_t(flip_i + 1) - flip_t(flip_i);
                    end
                    
                    % Load Stimulus                    
                    ev_e = new_pf.rec(i).ev_e{flip_ev_i(flip_i)};
                    stimulusParams = GratRevUtil.trigger2stimulusParams(GRW.pf, ev_e);
                    stimulus = GRW.stimulusLoader.getByStimulusParams(stimulusParams);
                    response = model.stimulate(stimulus);
                    
                    temporal_impulse_response = zeros(1000,1);
                    
%                     % Peak + Sustained Response
%                     transient_up = 0.2:0.2:6;
%                     transient_down = flip(transient_up);
%                     temporal_impulse_response(1:30) = transient_up;
%                     temporal_impulse_response(31:50) = transient_down(1:20);
%                     temporal_impulse_response(51:end) = transient_down(21);
%                     temporal_impulse_response = temporal_impulse_response / 4;

                    % Peak Only
                    PEAK_WIDTH = 40;
                    PEAK_HALF_WIDTH = round(PEAK_WIDTH/2);
                    transient_up = linspace(0,1,PEAK_HALF_WIDTH);
                    transient_down = flip(transient_up);
                    temporal_impulse_response(1:PEAK_HALF_WIDTH) = transient_up;
                    temporal_impulse_response(PEAK_HALF_WIDTH+1:2*PEAK_HALF_WIDTH) = transient_down;
                    %temporal_impulse_response(51:end) = transient_down(21);
                    
                    
                    transient_rate = 40;
                    
                    % Iterate ms by ms, generating spikes
                    for t=1:duration
                        % @TODO move poisson spike generation code into
                        % static method of ComplexCellModel
                        
                        % Add Noise to Firing Rate
                        firing_rate = (temporal_impulse_response(t) * (transient_rate + 0.001*response.firing_rate)) * (1 + 0.001*randn());
                        
                        % Generate Spike via Poisson Process
                        if rand() < (GRW.DT * firing_rate)
                            spike_time = flip_t(flip_i) + (t - 1) + GRW.LATENCY;
                            new_pf.rec(i).spike_times = [new_pf.rec(i).spike_times; spike_time];
                        end
                    end
                    
                    model.resetCache();
                end
            end
        end
            
        function new_pf = writeFromModel2(GRW, model, a, b)
        % WRITEFROMMODEL2 generates a new p2m file, NEW_PF, using the stimuli from
        % GRW.pf and a MODEL. The responses of the MODEL are passed through
        % a nonlinearity with terms A and B.
        %
        % NB: The model must have a method called stimulate(), which returns a
        % struct with a property called spike which is 0 or 1.
        
            new_pf = GRW.pf;
            
            % Outline:
            % 1. Iterate over each PF.rec
            % 2. Erase its spike_times
            % 3. Iterate over each frameflip event
            % 4. Load stimulus image
            % 5. Iterate ms by ms until next gap_frameflip evenet. If 
            %    spike occurs, add current time + LATENCY to spike_times
            
            for i=1:length(new_pf.rec)
                new_pf.rec(i).spike_times = [];
                
                [flip_ev_i, flip_t] = p2mFindEvents(new_pf,i,'FLIP');
                
                % Iterate over each flip
                for flip_i=1:length(flip_ev_i)
                    
                    % Calculate time interval between flips
                    if flip_i == length(flip_ev_i)
                        % if its the last flip, use default duration
                        % @todo calculate intelligently
                        duration = 100; 
                    else
                        duration = flip_t(flip_i + 1) - flip_t(flip_i);
                    end

                    % Load Stimulus                    
                    ev_e = new_pf.rec(i).ev_e{flip_ev_i(flip_i)};
                    stimulusParams = GratRevUtil.trigger2stimulusParams(ev_e);
                    stimulus = GRW.stimulusLoader.getByStimulusParams(stimulusParams);
                    response = model.stimulate(stimulus);
                    response_index = abs(response.SC1_response_sum) + abs(response.SC2_response_sum);
                    firing_rate = a*response_index^2 + b;
                    
                    % Iterate ms by ms, generating spikes
                    for t=1:duration
                        % @TODO move poisson spike generation code into
                        % static method of ComplexCellModel
                        
                        % Add Noise to Firing Rate
                        firing_rate = firing_rate * (1 + 0.1*randn());
                        
                        % Generate Spike via Poisson Process
                        if rand() < (GRW.DT * firing_rate)
                            spike_time = flip_t(flip_i) + (t - 1) + GRW.LATENCY;
                            new_pf.rec(i).spike_times = [new_pf.rec(i).spike_times; spike_time];
                        end
                    end
                    
                    model.resetCache();
                end
            end
        end
    end
end