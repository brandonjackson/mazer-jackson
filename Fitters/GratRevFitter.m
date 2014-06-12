classdef GratRevFitter < handle
    % GRATREVFITTER fits a model to a gratrev kernel
    % Each model has two custom methods in this class that is used to fit
    % its parameters to the gratrev kernel: one high-level function and
    % then an objective function (that is usuallly static) to use for
    % parameter optimization.
    %
    % NB: The model must have a method called stimulate(), which returns a
    % struct with a property called spike which is 0 or 1.

    
    properties
        stimulusLoader
        pf
        latency
        winsize
        tuning
        rast
        oriList
        sfList
        phaseList
        stimulusSize
        stimuli
        stimuliParams
        observedResponses
    end
    
	properties (Constant)
        LATENCY = 70;
    end
    
    methods
        function GRF = GratRevFitter(pf, latency, winsize)
        % GRATREVFITTER constructs an object that uses a MODEL to write
        % to PF, a p2m file.
        %
        % PARAMS
        %    pf
        %    model

            GRF.pf = pf;
            GRF.latency = latency;
            GRF.winsize=  winsize;
            
            % Get stimulus parameters
            GRF.tuning = pgratrev(GRF.pf, GRF.latency, GRF.winsize,'plot',0,'newfig',0);
            GRF.rast = prast(GRF.pf,'pre',GRF.latency,'post',GRF.latency+GRF.winsize);
            GRF.oriList = GRF.tuning{1,1}.x;
            GRF.sfList = GRF.tuning{2,1}.x;
            GRF.phaseList = GRF.tuning{3,1}.x;
            
            GRF.stimulusLoader = GratRevUtil.getStimulusLoader(GRF.pf);
            GRF.stimulusSize = GratRevUtil.getStimulusSize(GRF.pf);
            
            % convert gratrev kernel into ydata vector
            [GRF.stimuliParams, GRF.observedResponses] = GratRevUtil.getResponses(GRF.pf, GRF.latency, GRF.winsize);
            GRF.stimuli = cell(length(GRF.stimuliParams),1);
            for s=1:length(GRF.stimuliParams)
                
                % to prevent redundancy, the stimuli loader scans all
                % previously loaded stimuli and makes sure this one has not
                % been saved yet. If it has then it stores the index of the
                % first index instead of re-loading it.
                for d=1:(s-1)
                    if isequal(GRF.stimuliParams{s},GRF.stimuliParams{d})
                        GRF.stimuli{s} = d;
                        break;
                    end
                end
                if isempty(GRF.stimuli{s})
                    GRF.stimuli{s} = GRF.stimulusLoader.getByStimulusParams(GRF.stimuliParams{s});
                end
            end
        end
        
        function [ori, wavelength, sigma, a, b] = fitComplexCell(GRF, iterations)
        % FITCOMPLEXCELL fits a complex cell model by comparing the
        % observed and predicted responses.
            
            N_STARTS = iterations;
            
            options = optimset('Display','off','TolFun',1e-10);
            

                        
            % Assumption-Free Bounds
            lb = [0 5 2 -100 -100];
            ub = [180 150 100 100 100];
            
            % Opinionated Bounds
            % lb = [20 80 50 -1 0];
            % ub = [50 150 150 1 20];
            
            xs = zeros(N_STARTS,5);
            rs = zeros(N_STARTS,1);
            resnorms = zeros(N_STARTS,1);
            
            parfor i=1:N_STARTS
                % randomized starting point
                start = [randi([lb(1) ub(1)]) randi([lb(2) ub(2)]) randi([lb(3) ub(3)]) randi([lb(4) ub(4)]) randi([lb(5) ub(5)])];
                %start = [90 75 75 1 0]
                [x,resnorm] = lsqcurvefit(@GratRevFitter.fitComplexCellObjective,...
                    start,... % [ori wavelength sigma a b]
                    GRF.stimuli,... % stimuli vector used to generate xdata
                    GRF.observedResponses,...
                    lb,...
                    ub,...
                    options);
                xs(i,:) = x;
                resnorms(i) = resnorm;
                
                fittedData = GratRevFitter.fitComplexCellObjective(x,GRF.stimuli);
                r = corr(GRF.observedResponses,fittedData);
                rs(i) = r;
                
                fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',resnorm,r,x(1),x(2),x(3),x(4),x(5));
            end
            
            [~,minI] = min(resnorms);
            fprintf('and the winner is...\n');
            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',resnorms(minI),rs(minI),xs(minI,1),xs(minI,2),xs(minI,3),xs(minI,4),xs(minI,5));
            
            best_x = xs(minI,:);
            
            ori = best_x(1);
            wavelength = best_x(2);
            sigma = best_x(3);
            a = best_x(4);
            b = best_x(5);
        end

        function [ori, wavelength, sigma, a, b] = fitComplexCell2(GRF, iterations)
        % FITCOMPLEXCELL2 fits a complex cell model by comparing the
        % observed and predicted gratrev tuning curves
        %
        % @deprecated: use GRF.fitComplexCell() instead

            
            N_STARTS = iterations;
            
            options = optimset('Display','off','TolFun',1e-10);
            
            % convert gratrev kernel into ydata vector
            ydata = [GRF.tuning{1,1}.y; GRF.tuning{2,1}.y];
            
%             lb = [0 5 2 -100 -100];
%             ub = [180 150 150 100 100];
             lb = [20 80 50 -1 0];
             ub = [50 150 150 1 20];
            
            xs = zeros(N_STARTS,5);
            rs = zeros(N_STARTS,1);
            resnorms = zeros(N_STARTS,1);
            
            for i=1:N_STARTS
                % randomized starting point
                start = [randi([lb(1) ub(1)]) randi([lb(2) ub(2)]) randi([lb(3) ub(3)]) randi([lb(4) ub(4)]) randi([lb(5) ub(5)])];
                %start = [90 75 75 1 0]
                [x,resnorm] = lsqcurvefit(@GratRevFitter.fitComplexCellObjective2,...
                    start,... % [ori wavelength sigma a b]
                    zeros(length(ydata)),... % ignored by the objective function
                    ydata,...
                    lb,...
                    ub,...
                    options);
                xs(i,:) = x;
                resnorms(i) = resnorm;
                
                fittedData = GratRevFitter.fitComplexCellObjective2(x,zeros(length(ydata)));
                r = corr(ydata,fittedData);
                rs(i) = r;
                
                fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',resnorm,r,x(1),x(2),x(3),x(4),x(5));
            end
            
            [~,minI] = min(resnorms);
            fprintf('and the winner is...\n');
            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',resnorms(minI),rs(minI),xs(minI,1),xs(minI,2),xs(minI,3),xs(minI,4),xs(minI,5));
            
            best_x = xs(minI,:);
            
            ori = best_x(1);
            wavelength = best_x(2);
            sigma = best_x(3);
            a = best_x(4);
            b = best_x(5);
        end
        
        function [ori, wavelength, sigma, a, b] = fitComplexCellGlobal(GRF)
        % FITCOMPLEXCELL fits a complex cell model by comparing the
        % observed and predicted responses.
                        
            options = psoptimset('Display','iter','TolFun',1e-10,'MeshAccelerator','On','UseParallel','always','MaxIter',1500,'CompletePoll','on','MeshExpansion',4);%,'SearchMethod',{@searchlhs,1,25});
            
                        
            % Assumption-Free Bounds
            lb = [0 2 2 -1000 -1000];
            ub = [180 150 150 1000 1000];
            
            % Opinionated Bounds
            % lb = [20 80 50 -1 0];
            % ub = [50 150 150 1 20];


            % randomized starting point
            start = [30 randi([lb(2) ub(2)]) randi([lb(3) ub(3)]) 1 0];
            %start = [90 75 75 1 0]
            
            objective_corr = @(x)(1 - corr(GRF.observedResponses, GratRevFitter.fitComplexCellObjective(x,GRF.stimuli)));
            objective_lsq = @(x)(sum((GRF.observedResponses - GratRevFitter.fitComplexCellObjective(x,GRF.stimuli)).^2));

            problem = {};
            problem.x0 = start;
            problem.objective = objective_lsq;
            problem.lb = lb;
            problem.ub = ub;
            problem.options = options;

            [x,fval] = patternsearch(problem);

            fittedData = GratRevFitter.fitComplexCellObjective(x,GRF.stimuli);
            r = corr(GRF.observedResponses,fittedData);

            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',fval,r,x(1),x(2),x(3),x(4),x(5));

            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            a = x(4);
            b = x(5);
        end
        
        function [ori, wavelength, sigma, a, b, population, scores] = fitComplexCellGA(GRF)
        % FITCOMPLEXCELLGA fits a complex cell model by comparing the
        % observed and predicted responses.
                        
            options = gaoptimset('Display','diagnose','UseParallel','always','MutationFcn', {@mutationadaptfeasible, 0.4},'PopulationSize',400,'Generations',20,'StallGenLimit',100);%,'PlotInterval',5);%,'SearchMethod',{@searchlhs,1,25});
                
            % Assumption-Free Bounds
            lb = [0 2 2 -1000 -1000];
            ub = [180 150 150 1000 1000];
            
            % Opinionated Bounds
            % lb = [20 80 50 -1 0];
            % ub = [50 150 150 1 20];


            % randomized starting point
            start = [30 randi([lb(2) ub(2)]) randi([lb(3) ub(3)]) 1 0];
            %start = [90 75 75 1 0]
            
            fitness_corr = @(x)(1 - corr(GRF.observedResponses, GratRevFitter.fitComplexCellObjective(x,GRF.stimuli)));
            fitness_lsq = @(x)(sum((GRF.observedResponses - GratRevFitter.fitComplexCellObjective(x,GRF.stimuli)).^2));

            problem = {};
            %problem.x0 = start;
            %problem.objective = @(x)1 - corr(ydata, GratRevFitter.fitComplexCellObjective(x,stimuli));
            problem.fitnessfcn = fitness_lsq;
            problem.lb = lb;
            problem.ub = ub;
            problem.options = options;
            problem.solver = 'ga';
            problem.nvars = 5;

            [x,fval,~,~,population,scores] = ga(problem);

            fittedData = GratRevFitter.fitComplexCellObjective(x,GRF.stimuli);
            r = corr(GRF.observedResponses,fittedData);

            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',fval,r,x(1),x(2),x(3),x(4),x(5));

            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            a = x(4);
            b = x(5);
        end
        
        %function [ori, wavelength, sigma, a, b] = fitComplexCellSimpleGA(GRF)
        function results = fitComplexCellSimpleGA(GRF)
        % FITCOMPLEXCELLSIMPLEGA fits a complex cell model by fitting the gabor params, not the
        % nonlinearity
                        
            options = gaoptimset('Display','diagnose','UseParallel','always','MutationFcn', {@mutationadaptfeasible, 0.4},'PopulationSize',100,'Generations',20,'StallGenLimit',3);%,'PlotInterval',5);%,'SearchMethod',{@searchlhs,1,25});
                
            % Assumption-Free Bounds
            lb = [0 2 2];
            ub = [180 150 150];
            
            % Opinionated Bounds
            % lb = [20 80 50 -1 0];
            % ub = [50 150 150 1 20];


            % randomized starting point
            start = [30 randi([lb(2) ub(2)]) randi([lb(3) ub(3)])];
            %start = [90 75 75 1 0]
            
            % Generate subunit sum vector, then fits best nonlinearity to it
            fitness_nl = @(x)(GratRevFitter.fitNonlinearity(GratRevFitter.fitComplexCellObjectiveLinear(x,GRF.stimuli),GRF.observedResponses));

            problem = {};
            %problem.x0 = start;
            %problem.objective = @(x)1 - corr(ydata, GratRevFitter.fitComplexCellObjective(x,stimuli));
            problem.fitnessfcn = fitness_nl;
            problem.lb = lb;
            problem.ub = ub;
            problem.options = options;
            problem.solver = 'ga';
            problem.nvars = 3;

            [x,fval,~,~,population,scores] = ga(problem);

            subunitSums = GratRevFitter.fitComplexCellObjectiveLinear(x,GRF.stimuli);
            [sse,a,b] = GratRevFitter.fitNonlinearity(subunitSums,GRF.observedResponses);
            fittedData = a*subunitSums.^2 + b;
            r = corr(GRF.observedResponses,fittedData);

            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',fval,r,x(1),x(2),x(3),a,b);

            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            
            results.ori = ori;
            results.wavelength = wavelength;
            results.sigma = sigma;
            results.a = a;
            results.b = b;
            results.sse = sse;
            results.predictedResponses = fittedData;
            results.observedResponses = GRF.observedResponses;
            results.corr = r;
            %a = x(4);
            %b = x(5);
        end     
            
        function [ori, wavelength, sigma] = fitComplexCellGlobalLinear(GRF)
        % FITCOMPLEXCELL fits a complex cell model by comparing the
        % observed and predicted responses.
                        
            options = psoptimset('Display','iter','TolFun',1e-10,'MeshAccelerator','On','UseParallel','always','MaxIter',1500,'CompletePoll','on','MeshExpansion',4);%,'SearchMethod',{@searchlhs,1,25});
     
            % Assumption-Free Bounds
            lb = [0 2 2];
            ub = [180 150 150];
            
            % Opinionated Bounds
            % lb = [20 80 50 -1 0];
            % ub = [50 150 150 1 20];


            % randomized starting point
            start = [randi([lb(1) ub(1)]) randi([lb(2) ub(2)]) randi([lb(3) ub(3)])];
            %start = [90 75 75 1 0]

            problem = {};
            problem.x0 = start;
            %problem.objective = @(x)1 - corr(ydata, GratRevFitter.fitComplexCellObjective(x,stimuli));
            problem.objective = @(x)(sum((GRF.observedResponses - GratRevFitter.fitComplexCellObjectiveLinear(x,GRF.stimuli)).^2));
            problem.lb = lb;
            problem.ub = ub;
            problem.options = options;

            [x,fval] = patternsearch(problem);

            fittedData = GratRevFitter.fitComplexCellObjectiveLinear(x,GRF.stimuli);
            r = corr(GRF.observedResponses,fittedData);

            fprintf('f(x)=%.1f\tr=%.2f\t\tori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\n',fval,r,x(1),x(2),x(3));

            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
           % a = x(4);
           % b = x(5);
        end

        
    end
    
    methods (Static)
        
        function [sse, a, b] = fitNonlinearity(input,output)
            [f,gof] = fit(input,output,'a*x^2+b','StartPoint',[1 0]);
            sse = gof.sse;
            a = f.a;
            b = f.b;
        end
        
        function fitdata = fitComplexCellObjective(x, stimuli)
            % f = @(x,xdata)x(1)*xdata.^x(2)+x(3);
            
            %% Create a Complex Cell Model
            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            
            load CX_fit_args;
            load gratrev_stimulus_space;
            
            CX_fit_args.orientation = ori;
            CX_fit_args.gabor_params.wavelength = wavelength;
            CX_fit_args.gabor_params.sigma = sigma;
            
            CX = ComplexCellModel(CX_fit_args);
            
            %% Generate xdata by Stimulating Model with Stimuli
            xdata = zeros(length(stimuli),1);
            for i=1:length(stimuli)
                % If stimuli is a number, it is actually a reference to the
                % index where the stimulus is stored. This is part of a
                % caching system.
                if isscalar(stimuli{i})
                    xdata(i) = xdata(stimuli{i});
                else
                    response = CX.stimulate(stimuli{i});
                    xdata(i) = response.subunit_sum;
                end
            end
            
            %% Apply Nonlinearity to Generate Fitted Data
            a = x(4);
            b = x(5);
            fitdata = a*xdata.^2 + b;
            
            % fprintf('ori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\t\tNL=[%.6f,%.6f]\n',x(1),x(2),x(3),x(4),x(5));
        end

        function fitdata = fitComplexCellObjectiveLinear(x, stimuli)
        % FITCOMPLEXCELLOBJECTIVELINEAR returns the raw subunit sums
        % without applying a nonlinearity
        
            % f = @(x,xdata)x(1)*xdata.^x(2)+x(3);
            
            %% Create a Complex Cell Model
            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            
            load CX_fit_args;
            load gratrev_stimulus_space;
            
            CX_fit_args.orientation = ori;
            CX_fit_args.gabor_params.wavelength = wavelength;
            CX_fit_args.gabor_params.sigma = sigma;
            
            CX = ComplexCellModel(CX_fit_args);
            
            %% Generate xdata by Stimulating Model with Stimuli
            xdata = zeros(length(stimuli),1);
            for i=1:length(stimuli)
                if isscalar(stimuli{i})
                    xdata(i) = xdata(stimuli{i});
                else
                    response = CX.stimulate(stimuli{i});
                    xdata(i) = response.subunit_sum;
                end
            end
            
            fitdata = xdata;
            
          %  fprintf('ori=%.1f\tlambda=%.1fpx\tsigma=%.1fpx\n',x(1),x(2),x(3));
      end
      
      function fitdata = fitComplexCellObjective2(x, xdata)
            % f = @(x,xdata)x(1)*xdata.^x(2)+x(3);
            
            ori = x(1);
            wavelength = x(2);
            sigma = x(3);
            a = x(4);
            b = x(5);
            
            load CX_fit_args;
            load gratrev_stimulus_space;
            
            CX_fit_args.orientation = ori;
            CX_fit_args.gabor_params.wavelength = wavelength;
            CX_fit_args.gabor_params.sigma = sigma;
            
            CX = ComplexCellModel(CX_fit_args);
            
            oriResponses = zeros(length(sfList),length(oriList));
            sfResponses = zeros(length(oriList),length(sfList));
            
            for i=1:length(oriList)
                for j=1:length(sfList)
                    stimulusParams = {};
                    stimulusParams.stimulusSize = CX.stimulusSize(1);
                    stimulusParams.ori = oriList(i);
                    stimulusParams.sf = sfList(j);
                    stimulusParams.phase = 0;
                    stimulus = CX.stimulusLoader.getByStimulusParams(stimulusParams);
                    response = CX.stimulate(stimulus);
                    oriResponses(j,i) = response.subunit_sum;
                    sfResponses(i,j) = response.subunit_sum;
                end
            end
            
            oriMeans = mean(oriResponses,1);
            sfMeans = mean(sfResponses,1);
            
            xdata = [oriMeans';sfMeans']; % turn into a giant column
            
            fitdata = a*xdata.^2 + b;
        end
       
        
    end
end