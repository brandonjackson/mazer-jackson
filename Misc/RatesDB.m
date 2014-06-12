classdef RatesDB < handle
    %RATESDB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data % struct
    end
    
    properties (Constant)
        FILEDIR = '/lab/data/'
    end
    
    methods
        
        function RDB = RatesDB(varargin)
            p = inputParser;
            addParamValue(p,'experNames',0);
            addParamValue(p,'files',0);
            addParamValue(p,'data',0);
            parse(p,varargin{:});
            
            RDB.data = cell(1,1);
            if iscell(p.Results.experNames)
                RDB.build('experNames',p.Results.experNames);
            elseif iscell(p.Results.files)
                RDB.build('files',p.Results.files);
            elseif p.Results.data ~= 0
                RDB.data = p.Results.data;
            else
                RDB.build();
            end 
        end
        
        function [] = build(RDB, varargin)
        % builds/refreshes database
        
            p = inputParser;
            addRequired(p,'RDB');
            addParamValue(p,'experNames',0);
            addParamValue(p,'files',0);
            parse(p,RDB,varargin{:});
        
            if p.Results.experNames ~= 0
                experNames = p.Results.experNames;
            else
                if ~iscell(p.Results.files)
                    files = dir([RatesDB.FILEDIR '*gratrev*.p2m']);
                else
                    files = p.Results.files;
                end
                experNames = cell(length(files),1);
                for i=1:length(files)
                    [~,exper_name,exper_no] = fileparts(files(i).name);
                    experNames{i} = [exper_name exper_no];
                end
            end
            
            %experNames = experNames(1:10);
            
            for i=1:length(experNames)
                %[~,exper_name,exper_no] = fileparts(files(i).name);
                experName = experNames{i};
                if ~isempty(RDB.getByExperName(experName));
                    continue;
                end
                
                % Create row for new data files
                row = RDB.buildRow(experName);
                if ~isempty(row)
                    if RDB.isDataEmpty()
                        RDB.data{1} = row;
                    else
                        n_rows = length(RDB.data);
                        RDB.data{n_rows+1} = row;
                    end
                end
                fprintf('%s...done.\n',experName);
            end
        end
        
        function row = buildRow(RDB, experName)
            try
                pf = pffind(experName);
                raster = prast(pf);
                
                row = {};
                row.src = pf.src;
                row.experName = experName;
                row.taskName = PFUtil.taskName(pf);
                
                if strcmp(row.taskName,'gratrev')
                    row.isGrating = GratRevUtil.isGrating(pf);
                    row.isBar = GratRevUtil.isBar(pf);
                end
                
                row.reps = length(raster.triggers) / length(unique(raster.triggers));
                row.stimPeriod = round(1000/PFUtil.stimulusCarrierFrequency(pf));
                
                % PSTH
                [rates,times] = RasterUtil.rates(pf);
                row.psthRates = nanmean(rates,1);
                row.psthTimes = times;
                
                % Normalized PSTH
                N_PERIODS = 2;
                N_SAMPLES = 200;
                [rates_n, times_n] = RasterUtil.rates(pf,-1,row.stimPeriod*N_PERIODS+1); % pre of -1 necessary to get t(1)=0
                rates_mean_n = nanmean(rates_n,1);
                
                % Resample to have N_SAMPLES
                old_x = times_n;
                new_x = min(times_n):max(times_n)/(N_SAMPLES-1):max(times_n);
                row.psthNormalizedRates = interp1(old_x,rates_mean_n,new_x);
                row.psthNormalizedTimes = interp1(old_x,times_n,new_x);
                row.psthNormalizedPhase = linspace(0,2*pi*N_PERIODS,N_SAMPLES);
                
            catch err
                warning(getReport(err));
                row = [];
            end
        end
        
        function empty = isDataEmpty(RDB)
            empty = isempty(RDB.data{1});
        end
        
        function list = listPropertyValues(RDB,property)
            list = cell(length(RDB.data),1);
            n_hits = 0;
            for i=1:length(RDB.data)
                try
                    list{i} = RDB.data{i}.(property);
                    n_hits = n_hits + 1;
                catch err
                    % warning(getReport(err));
                end
            end
            if n_hits==0
                error(sprintf('Property `%s` not found',property));
            end
            if isnumeric(list{1}) && isnumeric(list{2})
                list = cell2mat(list);
            end
        end
        
        function row = getByExperName(RDB, experName)
            if RDB.isDataEmpty()
                row = [];
                return;
            end
            for i=1:length(RDB.data)
                if strcmp(RDB.data{i}.experName,'experName')==1
                    row = RDB.data{i};
                    return;
                end
            end
            row = []; % return empty if none found
        end
    end
    
    methods (Static)
        
        function results = filterByTaskName(data,taskName)
            if isempty(data{1})
                results = [];
                return;
            end
            n_matches = 0;
            results = cell(1,1);
            for i=1:length(data)
                if strcmp(data{i}.taskName,taskName)==1
                    n_matches = n_matches + 1;
                    results{n_matches} = data{i};
                end
            end
        end
        
        function results = filterByStimPeriod(data,stimPeriod)
            if isempty(data{1})
                results = [];
                return;
            end
            n_matches = 0;
            results = cell(1,1);
            for i=1:length(data)
                if (length(stimPeriod)==1 && data{i}.stimPeriod==stimPeriod)...
                   || (length(stimPeriod)==2 && data{i}.stimPeriod >= stimPeriod(1) && data{i}.stimPeriod <= stimPeriod(2))
                    n_matches = n_matches + 1;
                    results{n_matches} = data{i};
                end
            end
        end
        
    end
    
end

