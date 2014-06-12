classdef SuperWriter < handle
    % SUPERWRITER writes spike times to p2m files. This modifies the old
    % file, often using predicted spikes from custom models.
    
    properties
        stimulusLoader
        pf
        stimulusSize
    end
    
	properties (Constant)
        DT = 0.001;
    end
    
    methods
        function SW = SuperWriter()
        end
        
        function new_pf = writeFromModel(SW, model)
        % WRITEFROMMODEL generates a new p2m file, NEW_PF, using the stimuli from
        % GRW.pf and a MODEL
        %
        % NB: The model must have a method called stimulate(), which returns a
        % struct with a property called spike which is 0 or 1.
        end
        
        function [new_pf, ratio] = cliqueify(SW, interval, cardinality)
        % CLIQUEIFY eliminates all spikes whose inter-spike intervals are
        % greater than some INTERVAL (in ms), leaving behind only "cliques"
        % of spikes. Finally, the size of these cliques is examined, and
        % all cliques with less than CARDINALITY spikes are removed.
                
            if nargin < 3
                cardinality = 2;
            end
            new_pf = SW.pf;
            
            n_spikes = 0;
            n_spikes_in_cliques = 0;
        
            for i=1:length(SW.pf.rec)
                spikes = SW.pf.rec(i).spike_times;
                dts = diff(spikes);

                % find matching intervals
                dts_matching = dts <= interval;

                % recover indices of matching spikes from dt intervals (one on
                % either side of index of match, since each dt corresponds to
                % a pair of spikes)
                clique_left = [dts_matching; false];
                clique_right = [false; dts_matching];
                cliques_i = clique_left | clique_right;
                
                % Filter by cardinality of clique
                % Loops over each spike, and keeps track of how many
                % spikes are in each clique and then resets spikes that
                % don't meet cardinality requirements.
                if cardinality > 2
                    in_clique = 0;
                    clique_size = 0;
                    for j=1:length(spikes)
                        if ~in_clique
                            
                            % Entering new clique
                            if cliques_i(j)
                               % fprintf('j=%d entering clique, spike_time=%d\n',j,spikes(cliques_i(j)));
                                in_clique = 1;
                                clique_size = 1;
                                
                            % Still not in a clique
                            else
                                continue;
                            end
                        else
                            % Still in a clique
                            if cliques_i(j) && ((spikes(j) - spikes(j-1)) <= interval)
                                clique_size = clique_size + 1;
                                %fprintf('j=%d still in clique, new clique_size=%d\n',j,clique_size);
                               % fprintf('     isi=%d\n',(spikes(j) - spikes(j-1)));

                            % Exited a clique
                            else
                                % Clique not big enough, so ignore spikes
                                if clique_size < cardinality
                                    %fprintf('j=%d clique not big enough\n',j);
                                    cliques_i(j-clique_size:(j-1)) = false;
                                else
                                    % fprintf('j=%d clique exited, sufficiently large\n',j);
                                    % Sufficiently large clique found!
                                end
                                
                                if cliques_i(j)
                                    clique_size = 1;
                                    in_clique = 1;
                                else
                                    clique_size = 0;
                                	in_clique = 0;
                                end
                            end
                        end
                    end
                    
                    if in_clique && clique_size < cardinality
                        %fprintf('erasing tail end\n');
                    	cliques_i(length(spikes) - clique_size + 1:length(spikes)) = false;
                    end
                end
                
                spikes_in_cliques = spikes(cliques_i);
    
                n_spikes = n_spikes + length(spikes);
                n_spikes_in_cliques = n_spikes_in_cliques + length(spikes_in_cliques);

                new_pf.rec(i).spike_times = spikes_in_cliques;
            end
            
            ratio = n_spikes_in_cliques/n_spikes;
            fprintf('interval=%dms\t%0.0f%% of spikes preserved\n',interval,ratio*100);
            
        end
        
        function new_pf = discardRandomSpikes(SW, discardProbability)
        % DISCARDRANDOMSPIKES eliminates random spikes.
        % For each spike, there is a chance (with p=DISCARDPROBABILITY) 
        % that the spike will be discarded.
        %
        % The ratio of the number of spikes in NEW_PF to the number of 
        % spikes in the original pf approaches (1 - DISCARDPROBABILITY).
        %
        % USAGE
        %   SW.discardRandomSpikes(0.25); % discards 25% of spikes
        
            new_pf = SW.pf;
            
            n_spikes_pre = 0;
            n_spikes_post = 0;
        
            for i=1:length(SW.pf.rec)
                spikes = SW.pf.rec(i).spike_times;
                keepers = spikes(rand(length(spikes),1) >= discardProbability);
                
                n_spikes_pre = n_spikes_pre + length(spikes);
                n_spikes_post = n_spikes_post + length(keepers);

                new_pf.rec(i).spike_times = keepers;
            end
            
            % fprintf('%0.0f%% of spikes preserved\n',(n_spikes_post/n_spikes_pre)*100);
            
        end

        function new_pf = poissonSpikes(SW, winsize, homogeneous_rate)
        % POISSONSPIKES replaces the observed spikes with a new spike train
        % generated using an inhomegeneous Poisson process driven by the
        % estimated firing rate.
        %
        % If HOMOEGENEOUS_RATE is provided, the spikes are generated by a
        % homogeneous Poisson process with the provided rate (in spike/s)
        %
        % DEPRECATED
        
            if nargin < 2
                winsize = 20;
            end
        
            new_pf = SW.pf;
            
            % @todo add option for this
            if nargin >= 3
                new_pf.rec = [SW.pf.rec, SW.pf.rec, SW.pf.rec, SW.pf.rec, SW.pf.rec, SW.pf.rec, SW.pf.rec, SW.pf.rec];
            end
        
            for i=1:length(new_pf.rec)
                spikes = new_pf.rec(i).spike_times;
                new_spikes = [];
                [ts,rates] = SuperUtil.spikes2rates(spikes,'winsize',winsize);
                for j=1:length(ts)
                    
                    % Use homogeneous Poisson process if rate provided
                    if nargin < 3
                        p_spike = rates(j) * SuperWriter.DT;
                    else
                        p_spike = homogeneous_rate * SuperWriter.DT;
                    end
                    
                    if rand() <= p_spike
                        new_spikes = [new_spikes; ts(j)];
                    end
                end
                new_pf.rec(i).spike_times = new_spikes;
            end
        end
                
        function new_pf = inhomogeneousPoisson(SW, winsize, refractory_period, varargin)
        % POISSONSPIKES replaces the observed spikes with a new spike train
        % generated using an inhomegeneous Poisson process driven by the
        % estimated firing rate. Additionally, an absolute REFRACTORY_PERIOD is
        % enforced if provided.
        
            p = inputParser;
            addRequired(p,'winsize');
            addOptional(p,'refractory_period', 1);
            addParameter(p,'repetitions',1);
            parse(p,winsize, refractory_period,varargin{:});
        
            if nargin < 2
                winsize = 20;
            end
        
            new_pf = SW.pf;
            new_pf.rec = repmat(new_pf.rec,[1 p.Results.repetitions]);
                    
            for i=1:length(new_pf.rec)
                spikes = new_pf.rec(i).spike_times;
                [ts,rates] = SuperUtil.spikes2rates(spikes,'winsize',winsize);
                p_spikes = rates * SuperWriter.DT;
                rands = rand([length(p_spikes) 1]);
                above_threshold = rands <= p_spikes;
                
                % Look for spike in last REFRACTORY_PERIOD ms
                if nargin >= 3 && refractory_period > 1
                    spike_is = find(above_threshold==1);
                    for j=1:length(spike_is)
                        if above_threshold(spike_is(j))==1
                            above_threshold((spike_is(j)+1):(spike_is(j)+refractory_period-1)) = 0;
                        end
                    end
                end
                new_spikes = ts(above_threshold);
                new_pf.rec(i).spike_times = new_spikes;
            end
        end
        
        function new_pf = homogeneousPoisson(SW, homogeneous_rate, refractory_period, varargin)
        % HOMOGENEOUSPOISSON replaces the observed spikes with a new spike 
        % train generated using a homogeneous Poisson process driven by the
        % HOMOEGENEOUS_RATE. Additionally, an absolute REFRACTORY_PERIOD is
        % enforced.
        
            p = inputParser;
            addRequired(p,'homogeneous_rate');
            addOptional(p,'refractory_period', 1);
            addParameter(p,'repetitions',1);
            parse(p, homogeneous_rate, refractory_period, varargin{:});

                
            new_pf = SW.pf;
            new_pf.rec = repmat(new_pf.rec,[1 p.Results.repetitions]);
                    
            for i=1:length(new_pf.rec)
                spikes = new_pf.rec(i).spike_times;
                [ts,~] = SuperUtil.spikes2rates(spikes);
                ts = (min(ts)-400):(max(ts)+400); % prevent edge effects
                p_spikes = ones(size(ts)) * homogeneous_rate * SuperWriter.DT;
                rands = rand(size(p_spikes));
                above_threshold = rands <= p_spikes;
                
                % Look for spike in last REFRACTORY_PERIOD ms
                if nargin >= 3 && refractory_period > 1
                    spike_is = find(above_threshold==1);
                    for j=1:length(spike_is)
                        if above_threshold(spike_is(j))==1
                            above_threshold((spike_is(j)+1):(spike_is(j)+refractory_period-1)) = 0;
                        end
                    end
                end
                new_spikes = ts(above_threshold);
                new_pf.rec(i).spike_times = new_spikes;
            end
        end
     
    end
end