classdef GratRev < handle
    %GRATREV analyzes data from gratrev trials
    %
    % At the moment the primary job of this class is to extract a list of
    % oris, spatial frequencies and phases from a p2m file
    
    properties
        % p2m file
        pf
        
        % latency
        latency
        
        % winsize
        winsize
        
        exper
        
        plot_title
        
        params
        
        raster
        
        oris
        sfs
        phases
    end
    
    methods
        function GR = GratRev(pf, latency, winsize)
        % GRATREV constructs GratRev class
            
            GR.pf = getpf(pf);
            GR.latency = latency;
            GR.winsize = winsize;
                
            if isempty(strfind(GR.pf.rec(1).params.X_version_info{2}, 'gratrev'))
                error('not gratrev file');
            end
            
            % name of experiment (e.g. romeo0284.curvplay.005)
            [~,exper_name,exper_no] = fileparts(GR.pf.src);
            GR.exper = [exper_name exper_no];
            GR.plot_title = [GR.exper ' ' num2str(GR.latency) '-' num2str(GR.latency + GR.winsize) 'ms'];
            
            GR.params = GR.pf.rec(1).params;
            
            % Get raster
            GR.raster = prast(GR.pf);
            
            % Extract lists of oris, sfs and phases
            [GR.oris,GR.sfs,GR.phases] = GratRevUtil.getStimulusSpace(GR.pf);
        end
        
        function [] = plotBestStimuli(GR)
        % PLOTBESTSTIMULI plots the stimuli with highest evoked responses
        % Sorts stimuli based on their evoked response, then plots the 
        % kernel, details of the response distribution, and a grid of the
        % best stimuli.
            
            N_STIMS = 24; % stims to plot
            
            [resps,trigs] = RasterUtil.responses(GR.pf, GR.latency, GR.winsize);
            
            loader = GratRevUtil.getStimulusLoader(GR.pf);
            
            figure();
            
            % Plot Kernel
            subplot(5,6,1:2);
            grk2(GR.pf,GR.latency,GR.winsize);
            title('Kernel (sf vs. ori)');
            
            % Plot graph of Responses
            subplot(5,6,3:4);
            plot(1:N_STIMS,resps(1:N_STIMS),'-k.');
            hold on;
            plot(N_STIMS+1:length(resps),resps(N_STIMS+1:end),'-k');
            title('Firing Rates');
            
            % Plot Response Distribution
            subplot(5,6,5:6);
            hist(resps);
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',[0 0 0],'EdgeColor','w')
            title('Firing Rate Histogram');
            
            % Plot Best Stimuli
            for i=1:N_STIMS
                subplot(5,6,i+6);
                imshow(loader.getByTrigger(trigs{i}));
                axis square;
                title(sprintf('%.0f s/s',resps(i)));
            end
            
            boxtitle([PFUtil.experName(GR.pf) ' Best Stimuli']);
        end
    end
    
    methods(Static)
        
        function [oris, ori_responses] = getOriTuning(pf, latency, winsize)
            tunings = pgratrev(pf,latency,winsize,'plot',0);
            oris = tunings{1}.x;
            ori_responses = tunings{1}.y;
        end
        
    end
end
