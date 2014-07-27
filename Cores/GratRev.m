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
    end
    
    methods(Static)
        
        function [oris, ori_responses] = getOriTuning(pf, latency, winsize)
            tunings = pgratrev(pf,latency,winsize,'plot',0);
            oris = tunings{1}.x;
            ori_responses = tunings{1}.y;
        end
        
    end
end
