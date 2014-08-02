classdef STRF < handle
    % STRF finds spatio-temporal receptive fields
    
    properties
        pf
        validate
        theta_size
        stim_feature
        stim_freq
        
        best_lambda
        lambdas
        lambdas_scores
        lambdas_scores_err
        lambdas_scores_raw
        
        stims_cell
        rates_cell
        
        rates
        rates_validation
        stims
        stims_validation
        X
        X_normal
        y
        X_raw
        X_scale_mean
        X_scale_std
        
        BIN_WIDTH
        N_LAGS
        times
        
    end
    
    properties (Constant)


        N_EIGS = 10;
        N_JACKKNIFE_DRAWS = 10;
        TRAINING_DATA = 0.9
        
        N_NEGATIVE_LAGS = 0
    end
    
    methods
        function S = STRF(pf, varargin)
            p = inputParser;
            addRequired(p,'pf');
            addParameter(p,'feature','binary');
            addParameter(p,'auto_lambda',0);
            addParameter(p,'bin_width',0);
            addParameter(p,'n_lags',0);
            addParameter(p,'validate',0); % 0 or 1
            parse(p,pf,varargin{:});
            
%             % Low-Resolution
%             S.BIN_WIDTH = 16;
%             S.N_LAGS = 14;
% 
%             % Medium-Resolution
%             S.BIN_WIDTH = 12;
%             S.N_LAGS = 20;
% 
%             % Really-Good-Resolution
%             S.BIN_WIDTH = 10;
%             S.N_LAGS = 24;
            
            S.pf = pf;
            S.stim_freq = PFUtil.stimulusCarrierFrequency(S.pf);
            S.validate = p.Results.validate || (p.Results.auto_lambda==1);
            S.stim_feature = p.Results.feature;
            
            if p.Results.bin_width ~= 0
                S.BIN_WIDTH = p.Results.bin_width;
            end
            
            if p.Results.n_lags ~= 0
                S.N_LAGS = p.Results.n_lags;
            end
            
            S.setup();
            
            if p.Results.auto_lambda==1
                S.best_lambda = S.findBestLambda();
            end
        end
        
        function [] = setup(S)
            N_TRIALS = length(S.pf.rec);
            
            if isempty(S.BIN_WIDTH) || isempty(S.N_LAGS) 
                % High-Resolution
                S.BIN_WIDTH = 10;
                S.N_LAGS = 30;
            end
            
            S.times = -(S.N_NEGATIVE_LAGS*S.BIN_WIDTH):S.BIN_WIDTH:((S.N_LAGS-S.N_NEGATIVE_LAGS-1)*S.BIN_WIDTH);
            
            if S.validate == 1
                shuffled_recs = randperm(N_TRIALS,N_TRIALS);
                n_training_recs = round(N_TRIALS * STRF.TRAINING_DATA);
                training_rec_ix = shuffled_recs(1:n_training_recs);
                validation_rec_ix = shuffled_recs(n_training_recs+1:end);
            else
                training_rec_ix = 1:N_TRIALS;
            end
            
            if S.validate == 1
                % Gets full cell array of stim and rates on first run, and then
                % simply grabs slices and concatenates in the future, reducing
                % redundancy
                if isempty(S.stims_cell) || isempty(S.rates_cell)
                    S.stims_cell = PFUtil.stimulusMatrices(S.pf, 'feature', S.stim_feature,'padding',1);
                    S.rates_cell = PFUtil.rateVectors(S.pf,'padding',1);
                end

                S.rates = PFUtil.ratesCell2Vector(S.rates_cell(training_rec_ix),S.BIN_WIDTH);
                S.rates = S.rates - mean(S.rates);
                S.stims = PFUtil.stimulusCell2Matrix(S.stims_cell(training_rec_ix),S.BIN_WIDTH);
            
            
                S.rates_validation = PFUtil.ratesCell2Vector(S.rates_cell(validation_rec_ix),S.BIN_WIDTH);
                S.rates_validation = S.rates_validation - mean(S.rates_validation);
                S.stims_validation = PFUtil.stimulusCell2Matrix(S.stims_cell(validation_rec_ix),S.BIN_WIDTH);
            else
                % OLD WAY --- uses lots of redundant calls to
                % stimulusMatrices() when validate==1, but if only one draw
                % is required then this is much faster since it doesn't
                % keep stim_cells in memory (which can use over 1GB of
                % space)
                S.rates = PFUtil.rateVectors(S.pf,'concatenate',1,'bin_width',S.BIN_WIDTH,'rec_i',training_rec_ix);
                S.rates = S.rates - mean(S.rates);
                S.stims = PFUtil.stimulusMatrices(S.pf,'concatenate',1,'bin_width',S.BIN_WIDTH, 'feature', S.stim_feature, 'rec_i', training_rec_ix);
    %             
    %             if S.validate == 1
    %                 S.rates_validation = PFUtil.rateVectors(S.pf,'concatenate',1,'bin_width',S.BIN_WIDTH, 'rec_i', validation_rec_ix);
    %                 S.rates_validation = S.rates_validation - mean(S.rates_validation);
    %                 S.stims_validation = PFUtil.stimulusMatrices(S.pf,'concatenate',1,'bin_width',S.BIN_WIDTH, 'feature', S.stim_feature, 'rec_i', validation_rec_ix);
    %             end                
            end
                        
            S.theta_size = [size(S.stims,1) S.N_LAGS];
                        
            % If theta is huge, recompute stims and rates vectors with
            % bigger bins to prevent all hell from breaking loose!
            if S.theta_size(1)*S.theta_size(2) > 8000
                warning('Theta too large for regression, increasing bin size...');
                S.BIN_WIDTH = 12;
                S.N_LAGS = 20;
                S.setup();
                return;
            end
            
            % not sure if it's the best idea to do this up front, but it
            % does take a long time
            [S.X, S.y, S.X_raw, S.X_scale_mean, S.X_scale_std] = S.getXY(S.stims, S.rates);
            S.X_normal = S.X' * S.X;
        end
        
        function [X,y,X_raw,scale_mean,scale_std] = getXY(S, stims, rates)
            T = size(stims,2);
            stim_dims = size(stims,1);
            N_WEIGHTS = S.N_LAGS*stim_dims;
            
            %stims = sum(stims,1); % collapse to a single stimulus dimension
            
            regr_is = S.N_LAGS:(T - S.N_LAGS);
            X = zeros(length(regr_is), N_WEIGHTS);
            y = rates(regr_is);

            for i=1:length(regr_is)
                j = regr_is(i);
                X(i,:) = reshape(stims(:,(j - S.N_LAGS + 1 + S.N_NEGATIVE_LAGS):j+S.N_NEGATIVE_LAGS), [1 N_WEIGHTS]);
            end
            
            X_raw = sparse(X);
            
            % Demean and normalize to unit variance each column, required for ridge regression
            scale_mean = mean(X,1);
            X = bsxfun(@minus, X, scale_mean);
            
            % replace zeros with ones to prevent NaNs from appearing in X
            scale_std = std(X,1); 
            scale_std(scale_std==0) = 1;
            
            X = bsxfun(@rdivide, X, scale_std);
        end
        
        function [theta, theta_prediction, offset_prediction] = regress(S,varargin)
            p = inputParser;
            addRequired(p,'S');
            addParameter(p,'scrambled', 0);
            addParameter(p,'lambda', 0);
            parse(p,S,varargin{:});
            
            if p.Results.scrambled == 1
                y = S.y(randperm(length(S.y),length(S.y)));
            else
                y = S.y;
            end
            
            [theta, theta_prediction, offset_prediction] = S.ridge(S.X, y, p.Results.lambda);
            
        end
        
        function r_squared = crossValidate(S,varargin)
            p = inputParser;
            addRequired(p,'S');
            addParameter(p,'lambda', 0);
            parse(p,S,varargin{:});
            
            [~, theta_prediction, offset_prediction] = S.regress('lambda',p.Results.lambda);
            if S.validate
                [~, rates_actual, X_validation] = S.getXY(S.stims_validation, S.rates_validation);
                rates_predicted = (X_validation * theta_prediction) + offset_prediction;
            end
            r_squared = corr(rates_actual,rates_predicted)^2;
        end
        
        function best_lambda = findBestLambda(S,varargin)
            
            N_LAMBDAS = 10;
            N_DRAWS = 10;
            
            lambdas = logspace(-3,0,N_LAMBDAS);% 0:0.1:20;
            scores = zeros(N_DRAWS, N_LAMBDAS);
            fprintf('findBestLambda(): drawing ');
            for i=1:N_DRAWS
                fprintf('%d...',i);
                if i > 1 % re-initialize if not first trial
                    S.setup();
                end
                for j=1:N_LAMBDAS
                    scores(i,j) = S.crossValidate('lambda',lambdas(j));
                end
            end
            fprintf('\n');
            scores_mean = mean(scores,1);
            scores_err = std(scores,1);
            
            [~,max_i] = max(scores_mean);
            best_lambda = lambdas(max_i);
            
            S.lambdas = lambdas;
            S.lambdas_scores = scores_mean;
            S.lambdas_scores_err = scores_err;
            S.lambdas_scores_raw = scores;
        end
            
        function [] = plotLambdas(S)
            semilogx(S.lambdas,S.lambdas_scores - S.lambdas_scores_err,'-','Color',[0.8 0.8 0.8]);
            hold on;
            semilogx(S.lambdas,S.lambdas_scores,'.-','Color',[0 0.66 0]);
            semilogx(S.lambdas,S.lambdas_scores + S.lambdas_scores_err,'-','Color',[0.8 0.8 0.8]);
            
            % find max
            [~,max_i] = max(S.lambdas_scores);
            scatter(S.lambdas(max_i),S.lambdas_scores(max_i),100,[0 0.66 0],'p','MarkerFaceColor',[0 0.66 0]);
            
            xlim([min(S.lambdas),max(S.lambdas)]);

            title('Lambda Performance Analysis');
            ylabel('Prediction Score (r^2)');
            xlabel('lambda');
        end

        
%         function theta = temporalRegression(pf, varargin)
%             
%             p = inputParser;
%             addRequired(p,'pf');
%             addParameter(p,'temporal', 0);
%             addParameter(p,'separable',0);
%             addParameter(p,'temporal_rf', 0);
%             parse(p,pf,varargin{:});
% 
%             rates = PFUtil.rateVectors(pf,'concatenate',1,'bin_width',S.BIN_WIDTH);
%             rates = rates - mean(rates);
% 
%             stims = PFUtil.stimulusMatrices(pf,'concatenate',1,'bin_width',S.BIN_WIDTH);
%             
%             T = size(stims,2);
%             n_stims = size(stims,1);
%             theta = zeros(1);
%             
%             stims = sum(stims,1); % collapse to a single stimulus dimension
%             
%             regr_is = S.N_LAGS:(T - S.N_LAGS);
%             X = zeros(length(regr_is), S.N_LAGS);
%             y = rates(regr_is);
%             for i=1:length(regr_is)
%                 j = regr_is(i);
%                 X(i,:) = stims((j - S.N_LAGS + 1):j);
%             end
%             X = flip(X,2);
%             
%             theta = (X'*X)\(X'*y);
%             
%             if nargout < 1
%                 S.plotTemporalRF(pf, theta);
%             end
%         end

        function [eigenvalues_mean, eigenvalues_err] = scrambledEigs(S, varargin)
        % Not to be confused with scrambled eggs!    :-)
            p = inputParser;
            addRequired(p,'S');
            addParameter(p,'lambda',0);
            parse(p,S,varargin{:});          
            
            N_DRAWS = 15;
            
            eigenvalues_raw = zeros(N_DRAWS,STRF.N_EIGS);
            for i=1:N_DRAWS
                theta = S.regress('scrambled',1,'lambda',p.Results.lambda);
                [~,d] = eigs(cov(theta),STRF.N_EIGS);
                eigenvalues_unscaled = diag(d);
                eigenvalues_raw(i,:) = eigenvalues_unscaled / sum(eigenvalues_unscaled);
            end
            
            eigenvalues_mean = mean(eigenvalues_raw,1) * 100;
            eigenvalues_err = std(eigenvalues_raw,1) * 100;
        end
        
        function [theta_disp, theta_pred, offset_pred] = ridge(S, X, y, lambda_ratio)
        % THETA_DISP is a matrix with coefficient weights optimized for
        % display
        % THETA_PRED is a vector that can be used to generate predictions
        % OFFSET_PRED is a scalar that can be added to Y_hat to completely
        % reverse effects of scaling
            
            original_scale = 1;
            
            % We scale lambda based on the mean of the diagonals of B
            B_mean_diag = mean(diag(S.X_normal));
            lambda = B_mean_diag * lambda_ratio;
            tikhonov_matrix = (lambda * eye(size(S.X_normal)));
            
            % Find the beta weights!
            theta_v = (S.X_normal + tikhonov_matrix) \ (X' * y);
            
            if original_scale == 1
                %theta_v_unscaled = theta_v;
                theta_pred = theta_v ./ S.X_scale_std';
                offset_pred = mean(y) - (S.X_scale_mean * theta_pred);
                %theta_v = theta_v_rescaled;
            end
            
%             y_hat_fake = X * theta_v_unscaled;
%             y_hat_real = (S.X_raw * theta_v) + offset;
%             
%             corr(y_hat_fake, y_hat_real)
            
            theta = reshape(theta_v, S.theta_size);
            theta_disp = flip(theta,2);
        end
        
        function [] = printBatchLambdas(S)
            
            N_LAMBDAS = 10;
            LAMBDAS = logspace(-1.5,1.5,N_LAMBDAS);
            
            [~,exper_name,exper_no] = fileparts(S.pf.src);
            
            for i=1:N_LAMBDAS
                lambda = LAMBDAS(i);
                filename = ['/lab/results/batch_pstrf/' exper_name exper_no '_pstrf_lambda_' sprintf('%.2f',lambda) '.png'];
                S.plotTheta(lambda, [' STRF (\lambda = ' sprintf('%.2f',lambda) ', feature=binary )']);
                sdf('printable');
                orient portrait;
                print('-dpng',filename);
                close all;
            end
            
        end
        
                
        function [theta_final, V, eigenvalues,sig_eigs] = plotTheta(S, lambda, detailsStr)
            
            if nargin < 3
                detailsStr = '';
            end
            
            pf = S.pf;
            
            % Disable validation, use full dataset
            old_validate = S.validate;
            if S.validate==1
                S.validate = 0;
                S.setup();
            end
            
            [theta, theta_prediction,offset_prediction] = S.regress('lambda',lambda);
            
%             theta_final = reshape(theta_prediction, S.theta_size);
%             theta_final = flip(theta_final,2);
            
            theta_final = theta;
            
            % PCA
            [V,eigenvalues] = S.pca(theta_final);
            
            % Sort theta by strength of projection onto first eigen vector
            projections = zeros(size(theta_final,1),1);
            for i=1:size(theta_final,1)
                projections(i) = theta_final(i,:) * V(:,1);
            end
            [~,theta_i] = sort(abs(projections),'descend');
            theta_final = theta_final(theta_i,:);
            
            % Scrambled PCA (to determine whether eigenvalues are
            % significant)
            [scrambled_eigs, scrambled_eigs_err] = S.scrambledEigs('lambda',lambda);
            
            sig_eigs = zeros(STRF.N_EIGS,2);
            
            for i=1:STRF.N_EIGS
                %[sig_eigs(i,1),sig_eigs(i,2)] = ztest(eigenvalues(i)*100,scrambled_eigs(i),scrambled_eigs_err(i),'Alpha',0.05,'Tail','right');
                [sig_eigs(i,1),sig_eigs(i,2)] = ztest(eigenvalues(i)*100,scrambled_eigs(i),scrambled_eigs_err(i),0.01,'right');
            end
            
            % count # sig eigen vectors, stopping once a non-significant is
            % encountered
            n_sig_eigs = 0;
            sig_eig_sum = zeros(1,size(theta_final,2));
            for i=1:STRF.N_EIGS
                if sig_eigs(i,1)==1
                    n_sig_eigs = n_sig_eigs + 1;
                    sig_eig_sum = sig_eig_sum + eigenvalues(i)*V(:,i)';
                else
                    break;
                end
            end
                
            
            raster = prast(S.pf);
            n_frames = length(raster.triggers);
            n_reps = n_frames / length(unique(raster.triggers));

            figure();

            subplot(3,3,1);
            S.plotTemporalRF(pf, V(:,1));
            title('Eigen Vector #1');
            xlabel('');
            
            subplot(3,3,2);
            S.plotTemporalRF(pf, V(:,2));
            title('Eigen Vector #2');
            xlabel('');
            
            subplot(3,3,4);
            errorbar(scrambled_eigs,scrambled_eigs_err,'Color',[0.8 0.8 0.8]);
            hold on;
            plot(eigenvalues*100,'-k.');
            
            
            star_height = max(eigenvalues*100)*1.2;
            star_size = 36;
            for i=1:STRF.N_EIGS
                if sig_eigs(i,1)==1
                    if sig_eigs(i,2) <= 0.001
                        scatter(i-0.2,star_height,star_size,'k','*');
                        scatter(i,star_height,star_size,'k','*');
                        scatter(i+0.2,star_height,star_size,'k','*');
                    elseif sig_eigs(i,2) <= 0.01
                        scatter(i-0.1,star_height,star_size,'k','*');
                        scatter(i+0.1,star_height,star_size,'k','*');
                    else
                        scatter(i,star_height,star_size,'k','*');
                    end
                else
                    break;
                end
            end
            
            xlim([0 length(eigenvalues)]);
            ylim([0 max(star_height*1.1,33)]);
            xlabel('n');
            ylabel('% variance explained');
            legend('null','actual','Location','NorthEast');
            title('Eigenvalues');

            subplot(3,3,5);
            if ~isempty(S.best_lambda)
                S.plotLambdas();
            else
                if n_sig_eigs > 0
                    S.plotTemporalRF(S.pf, sig_eig_sum);
                    title('Scaled Sum of Significant Eigen Vectors');
                else
                    title('No Significant Eigen Vectors');
                end
            end
            
            subplot(3,3,3);
            S.plotTemporalRF(pf, mean(theta_final,1));
            title('Mean Row Vector');
            
            subplot(3,3,7);
            PFUtil.plotCarrierFrequency(S.pf,0);
            xlim([min(S.times) max(S.times)]);
            set(gca,'XGrid','off');
%             [stim_freq, stim_dur, stim_gapdur] = PFUtil.stimulusCarrierFrequency(S.pf);
%             
%             psth = phist(pf,'pre',0,'post',S.BIN_WIDTH * S.N_LAGS,'suptitle',0);
%             ts = psth(:,1);
%             rates = psth(:,2);
%             variance = psth(:,3);
%             eshade(ts,rates,variance);
%             hold on;
%             plot(ts,rates,'-r');
%             
%             y_max = 1.1*(max(rates)+max(variance));
% 
%             
%             if stim_gapdur > 0
%                 line([stim_dur,stim_dur],[0,y_max],'Color',[0.5 0.5 0.5]);
%                 line([stim_dur+stim_gapdur,stim_dur+stim_gapdur],[0,y_max],'Color','k');
%             else
%                 n_stim_onsets = floor((S.BIN_WIDTH * S.N_LAGS)/stim_dur);
%                 for i=1:n_stim_onsets
%                     line([i*stim_dur,i*stim_dur],[0,y_max],'Color','k');
%                 end
%             end
%             xlim([0 S.BIN_WIDTH * S.N_LAGS]);
%             ylim([0 y_max]);           
%             title(sprintf('PSTH  (%d frames, %.1f reps)',n_frames,n_reps));
%             ylabel('s/s');
%             xlabel('Lag (ms)');
            
            subplot(3,3,8);
            if n_sig_eigs > 0
                step = zeros(((S.BIN_WIDTH * S.N_LAGS) / S.BIN_WIDTH)+10, 1);
                step(11:end) = 1;
                step_resp = conv(step,sig_eig_sum);
                onset = find(step_resp,1,'first');
                plot(step_resp(onset+1:onset+S.N_LAGS),'-r');
                hold on;
                line([0 S.N_LAGS],[0 0],'Color','k');
                xlim([0 S.N_LAGS]);
                ylim([0 1.3*max(step_resp)]);
                title('Significant Eigen Vectors Step Response');
                xticks = 5:5:S.N_LAGS;
                legendCell = cellstr(num2str((xticks*S.BIN_WIDTH)'));
                set(gca,'XTick',xticks);
                set(gca,'XTickLabel',legendCell);
                xlabel('Lag (ms)');
            else
                title('No Significant Eigen Vectors');
            end
            
            subplot(3,3,[6,9]);
            imagesc(theta_final);
            colormap(hotcold(100));
            %colorbar();
            xlabel('Lag (ms)');
            ylabel('Stimulus Features');
            set(gca,'YTick',[]);
            xticks = 5:10:S.N_LAGS;
            legendCell = cellstr(num2str(((xticks*S.BIN_WIDTH)-(S.N_NEGATIVE_LAGS*S.BIN_WIDTH))'));
            set(gca,'XTick',xticks);
            set(gca,'XTickLabel',legendCell);
            zmax = max(max(theta_final(:)),abs(min(theta_final(:))));
            caxis([-zmax zmax]);
            title('Beta Matrix');
            
            experName = PFUtil.experName(S.pf);
            
            if strcmp(PFUtil.taskName(S.pf),'gratrev')
                if GratRevUtil.isBar(S.pf)
                    experName = [experName ' (Bars)'];
                else
                    experName = [experName ' (Gratings)'];
                end
            end
            
            titleStr = sprintf('%s STRF (lambda = %.2f, feature=%s, %d %dms bins)',...
            experName,...
            lambda,...
            S.stim_feature,...
            S.N_LAGS,...
            S.BIN_WIDTH);
            boxtitle(titleStr);
            
            %boxtitle([PFUtil.experName(pf) ' ' detailsStr]);
            
            S.validate = old_validate;
        end
        
        function [theta_final, theta_err, impulse] = plotThetaDelta(S, lambda, detailsStr)
            
            if nargin < 3
                detailsStr = '';
            end
            
            pf = S.pf;
            
            % Disable validation, use full dataset
            old_validate = S.validate;
            if S.validate==1
                S.validate = 0;
                S.setup();
            end
            
            [theta,theta_err] = PFUtil.jackknife(S.pf,@STRF.jackknifeWrapper,10,lambda,'bin_width',S.BIN_WIDTH,'n_lags',S.N_LAGS);
            impulse = [0;diff(theta)];
            %[theta, theta_prediction,offset_prediction] = S.regress('lambda',lambda);
            
%             theta_final = reshape(theta_prediction, S.theta_size);
%             theta_final = flip(theta_final,2);
            
            theta_final = theta;
            
            % PCA
%             [V,eigenvalues] = S.pca(theta_final);
%             
%             % Sort theta by strength of projection onto first eigen vector
%             projections = zeros(size(theta_final,1),1);
%             for i=1:size(theta_final,1)
%                 projections(i) = theta_final(i,:) * V(:,1);
%             end
%             [~,theta_i] = sort(abs(projections),'descend');
%             theta_final = theta_final(theta_i,:);
%             
%             % Scrambled PCA (to determine whether eigenvalues are
%             % significant)
%             [scrambled_eigs, scrambled_eigs_err] = S.scrambledEigs('lambda',lambda);
%             
%             sig_eigs = zeros(STRF.N_EIGS,2);
%             
%             for i=1:STRF.N_EIGS
%                 %[sig_eigs(i,1),sig_eigs(i,2)] = ztest(eigenvalues(i)*100,scrambled_eigs(i),scrambled_eigs_err(i),'Alpha',0.05,'Tail','right');
%                 [sig_eigs(i,1),sig_eigs(i,2)] = ztest(eigenvalues(i)*100,scrambled_eigs(i),scrambled_eigs_err(i),0.01,'right');
%             end
            
            % count # sig eigen vectors, stopping once a non-significant is
            % encountered
            n_sig_eigs = 0;
%             sig_eig_sum = zeros(1,size(theta_final,2));
%             for i=1:STRF.N_EIGS
%                 if sig_eigs(i,1)==1
%                     n_sig_eigs = n_sig_eigs + 1;
%                     sig_eig_sum = sig_eig_sum + eigenvalues(i)*V(:,i)';
%                 else
%                     break;
%                 end
%             end
                
            
            raster = prast(S.pf);
            n_frames = length(raster.triggers);
            n_reps = n_frames / length(unique(raster.triggers));

            figure();
            if ~isempty(S.best_lambda)
                N_ROWS = 3;
            else
                N_ROWS = 2;
            end
            
            
            subplot(N_ROWS,2,1);
            S.plotTemporalRF(pf, theta_final, theta_err);
            title('Beta Weights + Jackknife stderr');
            
            subplot(N_ROWS,2,2);
            S.plotTemporalRF(pf, impulse);
            title('Diff of Betas / Temporal Impulse Response');

            subplot(N_ROWS,2,3);
            PFUtil.plotCarrierFrequency(S.pf,0);
            xlim([min(S.times) max(S.times)]);
            set(gca,'XGrid','off');
            
            subplot(N_ROWS,2,4);
            step = zeros(S.N_LAGS+10, 1);
            step(11:end) = 1;
            step_resp = conv(step,impulse);
            onset = find(step_resp,1,'first');
            plot(step_resp(onset+1:onset+S.N_LAGS),'-r');
            hold on;
            
%             line([0 max(S.times)],[0 0],'Color','k');
%             xlim([min(S.times) max(S.times)]);
            
            line([0 S.N_LAGS],[0 0],'Color','k');
            xlim([0 S.N_LAGS]);
            
            ylim([0 1.3*max(step_resp)]);
            title('Impulse Step Response');
            
            xticks = 5:10:S.N_LAGS;
            legendCell = cellstr(num2str(((xticks*S.BIN_WIDTH)-(S.N_NEGATIVE_LAGS*S.BIN_WIDTH))'));
            
%             xticks = 5:5:S.N_LAGS;
%             legendCell = cellstr(num2str((xticks*S.BIN_WIDTH)'));
            
            set(gca,'XTick',xticks);
            set(gca,'XTickLabel',legendCell);
            xlabel('Lag (ms)');
            
            
            if ~isempty(S.best_lambda)
                subplot(N_ROWS,2,5);
                S.plotLambdas();
            end
            
            experName = PFUtil.experName(S.pf);
            
            if strcmp(PFUtil.taskName(S.pf),'gratrev')
                if GratRevUtil.isBar(S.pf)
                    experName = [experName ' (Bars)'];
                else
                    experName = [experName ' (Gratings)'];
                end
            end
            titleStr = sprintf('%s STRF (lambda = %.2f, feature=%s, %d %dms bins)',...
                experName,...
                lambda,...
                S.stim_feature,...
                S.N_LAGS,...
                S.BIN_WIDTH);
            boxtitle(titleStr);
            
            S.validate = old_validate;
        end
        
         function [] = plotImageKernel(S,theta)
            n_pixels = size(theta,1);
            n_lags = size(theta,2);
            subplot_width = ceil(sqrt(n_lags));
            if (subplot_width * (subplot_width - 1)) >= n_lags
                subplot_height = subplot_width - 1;
            else
                subplot_height = subplot_width;
            end
            color_scale = [-max(abs(theta(:))) max(abs(theta(:)))];
            
            for i=1:n_lags
                subplot(subplot_height, subplot_width, i);
                img = reshape(theta(:,i),[sqrt(n_pixels) sqrt(n_pixels)]);
                imagesc(img);
                caxis(color_scale);
                set(gca,'YTick',[]);
                set(gca,'XTick',[]);
                axis square;
                title(['Lag = ' num2str((i-1)*S.BIN_WIDTH) 'ms']);
            end
         end
        
         function [] = plotTemporalRF(S, pf, rf, err)
            xs = S.times;
            
            if nargin == 4
                errorbar(xs,rf,err,'-b');
            else
                plot(xs, rf, '-b.');
            end
            y_max = max(abs(rf))*1.2;
            hold on;
            line([min(xs) max(xs)],[0 0],'Color','k');
            
            for i=1:5
                onset_t = (i-1)*1000*(1/S.stim_freq);
                line([onset_t,onset_t],[-y_max y_max],'Color',[0.5 0.5 0.5]);
            end
            title([PFUtil.experName(pf) ' Temporal Impulse Response']);
            xlabel('Lag (ms)');
            xlim([min(xs) max(xs)]);
            if y_max==0
                ylim([-0.5 0.5]);
            else
                ylim([-y_max y_max]);
            end
         end
                 
        function [eigen_vectors,eigenvalues] = pca(S,theta)
            
            [eigen_vectors,d] = eigs(cov(theta), STRF.N_EIGS);
            eigenvalues_raw = diag(d);
            eigenvalues = eigenvalues_raw / sum(eigenvalues_raw);
            
            % The sign of the eigen vectors is arbitrary. To ensure
            % consistency, though, we will flip the sign so that the area
            % beneath each curve (between 60 and 100) is greater than zero.
            min_i = floor(60 / S.BIN_WIDTH);
            max_i = floor(100 / S.BIN_WIDTH);
            for i=1:STRF.N_EIGS
                area = trapz(eigen_vectors(min_i:max_i,i));
                if area < 0
                    eigen_vectors(:,i) = eigen_vectors(:,i)*-1;
                end
            end
        end
        
%         function [theta_mean,theta_err] = jackknife(S,lambda,varargin)
% 
%             pf = S.pf;
%             
%             N_RECS = length(pf.rec);
%             
%             shuffled_recs = randperm(N_RECS);
%             recs_per_draw = floor(N_RECS / n_draws);
%             
%             results = cell(n_draws,1);
%             max_rec_i = 0;
%             for i=1:n_draws
%                 rand_pf = pf;
%                 if i~=n_draws
%                     shuffled_recs_i = (max_rec_i+1):(max_rec_i + recs_per_draw);
%                 else
%                     shuffled_recs_i = (max_rec_i+1):length(shuffled_recs);
%                 end
%                 excluded_recs = shuffled_recs(shuffled_recs_i);
%                 remaining_recs = setdiff(1:N_RECS, excluded_recs);
%                 rand_pf.rec = rand_pf.rec(remaining_recs);
%                 results{i} = f(rand_pf,varargin{:});
%                 max_rec_i = max(shuffled_recs_i);
%             end
%             
%             result_size = size(results{1});
%             results_matrix = zeros([n_draws, result_size]);
%             
%             for i=1:n_draws
%                 if min(result_size)==1
%                     results_matrix(i,:) = results{i};
%                 else
%                     results_matrix(i,:,:) = results{i};
%                 end
%             end
%             
%             result_mean = squeeze(nanmean(results_matrix,1));
%             result_std = squeeze(nanstd(results_matrix,1));
%         end
        
    end
    
    methods (Static)

        function [] = meanSTRF(stimulus,lambda,n_clusters)

            path = ['/lab/results/batch_pstrfs/' stimulus '/'];
            files_query = sprintf('%s*%.2f*.mat',path,lambda);
            files = dir(files_query);

            if strcmp(stimulus,'angleplay')
                N_LAGS = 24;
            else
                N_LAGS = 25;
            end

            vectors = zeros(length(files),N_LAGS);
            lags = [];
            for i=1:length(files)
                filename = [path files(i).name];
                load(filename);

                P = pstrf_result;
                lags = P.lags(1:N_LAGS);
                vector_sig_acc = zeros(N_LAGS,1);
                n_sigs = 0;
                for j=1:length(P.eigenvalues_significance)
                    if P.eigenvalues_significance(j,1) == 1
                       vector_sig_acc = vector_sig_acc + (P.eigen_vectors(:,j)*P.eigenvalues(j));
                       n_sigs = n_sigs + 1;
                    end
                end

                weighted_eigen_vector = vector_sig_acc / n_sigs;
                first_eigen_vector = P.eigen_vectors(:,1);
                mean_beta = P.theta_mean;

                % Default Scaling (i.e. NONE!)
                vectors(i,:) = mean_beta;

                % Scale Max to 1
                % vectors(i,:) = vectors(i,:) / max(abs(vectors(i,:)));

                % Scale to Unit Variance
                vectors(i,:) = vectors(i,:) / std(vectors(i,:));
            end

            % get rid of nan vectors
            vectors = vectors(sum(isnan(vectors),2)==0,:);

            % K-means clustering
            N_SUBPLOTS = n_clusters + 2;
            N_COLS = 2;
            N_ROWS = ceil(N_SUBPLOTS/N_COLS);
            [cluster_ix,clusters] = kmeans(vectors, n_clusters);
            color_palette = cell(4,1);
            color_palette{1} = 'r';
            color_palette{2} = 'b';
            color_palette{3} = 'g';
            color_palette{4} = 'c';

            % MDS
            dissimilarities = pdist(vectors);
            [MDS_Y,stress] = mdscale(dissimilarities,2,'criterion','metricstress');

            figure();

            subplot(N_ROWS,N_COLS,1);
            line([min(lags) max(lags)],[0 0],'Color','k');
            hold on;
            errorbar(lags,mean(vectors,1),std(vectors,1),'Color',[0.9 0.9 0.9]);
            plot(lags,mean(vectors,1),'LineWidth',1.5,'Color','k');
            xlim([min(lags) max(lags)]);
            title('Mean Temporal Impulse Response');
            xlabel('lag (ms)');

            subplot(N_ROWS,N_COLS,2);

            for i=1:n_clusters
                scatter(MDS_Y((cluster_ix==i),1),MDS_Y((cluster_ix==i),2),36,color_palette{i});
                hold on;
            end
            title('Multidimensional Scaling');

            max_cluster_y = max(clusters(:))*1.2;
            for i=1:n_clusters
                subplot(N_ROWS,N_COLS,i+2);
                line([min(lags) max(lags)],[0 0],'Color','k');
                hold on;
                plot(lags,clusters(i,:),'LineWidth',1.5,'Color',color_palette{i});
                xlim([min(lags) max(lags)]);
                ylim([-max_cluster_y max_cluster_y]);
                xlabel('lag (ms)');
                title(sprintf('K-Means Centroid #%d (n=%d)',i,sum(cluster_ix==i)));
            end
        end
        
        function [] = meanDeltaSTRF(path, stimulus,lambda,n_clusters,min_score,plot_kernels)
            
            if nargin < 5
                min_score = 0.1;
            end
            
            if nargin < 6
                plot_kernels = 0;
            end
            
           % path = ['/lab/results/batch_pstrf_deltas/' stimulus '/'];
            if strcmp(lambda,'auto')
                files_query = sprintf('%s*auto*.mat',path);
            else
                files_query = sprintf('%s*%.2f*.mat',path,lambda);
            end
            files = dir(files_query);

            if strcmp(stimulus,'angleplay') || strcmp(stimulus,'gridcurv') || strcmp(stimulus,'reel')
                N_LAGS = 24;
            else
                N_LAGS = 15;
            end

            vectors = zeros(length(files),N_LAGS);
            lags = [];
            
            if plot_kernels==1
                figure();
                kernel_plot_i = 1;
            end
            
            for i=1:length(files)
                filename = [path files(i).name];
                load(filename);

                P = pstrf_result;
                lags = P.lags(1:N_LAGS);
                
                if max(P.lambdas_scores) < min_score
                    continue;
                end
                

                % Potential vectors to use for meta-analysis
                mean_beta = P.theta;
                impulse = P.theta_diff;
                impulse(1:2) = 0;
                impulse(end) = 0;

                % Default Scaling (i.e. NONE!)
                vectors(i,:) = impulse;

                % Scale Max to 1
                % vectors(i,:) = vectors(i,:) / max(abs(vectors(i,:)));

                % Scale to Unit Variance
                vectors(i,:) = vectors(i,:) / std(vectors(i,:));
                
                if plot_kernels && kernel_plot_i <= 35
                    subplot(5,7,kernel_plot_i);
                    plot(impulse,'-r');
                    hold on;
                    plot(mean_beta,'-b');
                    title(P.src);
                    kernel_plot_i = kernel_plot_i + 1;
                end
            end

            % get rid of nan vectors
            n_vectors_original = size(vectors,1);
            vectors = vectors(sum(isnan(vectors),2)==0,:);
            vectors = vectors(abs(sum(vectors,2)) > 0, :);
            n_vectors_final = size(vectors,1);
            fprintf('STRF.meanDeltaSTRF(): %d / %d discarded due to poor prediction scores\n',...
                (n_vectors_original - n_vectors_final), n_vectors_original);

            % K-means clustering
            N_SUBPLOTS = n_clusters + 2;
            N_COLS = 2;
            N_ROWS = ceil(N_SUBPLOTS/N_COLS);
            [cluster_ix,clusters] = kmeans(vectors, n_clusters);
            color_palette = cell(6,1);
            color_palette{1} = 'r';
            color_palette{2} = 'b';
            color_palette{3} = 'g';
            color_palette{4} = 'c';
            color_palette{5} = 'm';
            color_palette{6} = [0.5 0.5 0.5];

%             %MDS
%             dissimilarities = pdist(vectors);
%             MDS_Y = mdscale(dissimilarities,2,'criterion','metricstress');

            figure();

            subplot(N_ROWS,N_COLS,1);
            line([0 max(lags)],[0 0],'Color','k');
            hold on;
            errorbar(lags,mean(vectors,1),std(vectors,1),'Color',[0.9 0.9 0.9]);
            plot(lags,mean(vectors,1),'LineWidth',1.5,'Color','k');
            xlim([0 max(lags)]);
            title('Mean Temporal Impulse Response');
            xlabel('lag (ms)');

%             subplot(N_ROWS,N_COLS,2);
% 
%             for i=1:n_clusters
%                 scatter(MDS_Y((cluster_ix==i),1),MDS_Y((cluster_ix==i),2),36,color_palette{i});
%                 hold on;
%             end
%             title('Multidimensional Scaling');

            max_cluster_y = max(clusters(:))*1.2;
            for i=1:n_clusters
                subplot(N_ROWS,N_COLS,i+2);
                line([0 max(lags)],[0 0],'Color','k');
                hold on;
                plot(lags,clusters(i,:),'LineWidth',1.5,'Color',color_palette{i});
                xlim([0 max(lags)]);
                ylim([-max_cluster_y max_cluster_y]);
                xlabel('lag (ms)');
                title(sprintf('K-Means Centroid #%d (n=%d)',i,sum(cluster_ix==i)));
            end
            
            boxtitle(sprintf('%c%s Mean Temporal IRs (min prediction = %.2f)',...
                upper(stimulus(1)),...
                stimulus(2:end),...
                min_score));
           
        end
        
        function theta = jackknifeWrapper(pf,lambda,varargin)
            S = STRF(pf,'feature','delta',varargin{:});
            theta = S.regress('lambda',lambda);
        end
        
        function [] = batch_pstrf_delta(file_list, lambda, filedir)
        % BATCH_PSTRF_DELTA runs and saves pstrf_delta results for many
        % cells
        %
        % INPUTS
        %   file_list - cell array with file names (or any thing that can
        %               be fed into dbfind)
        %   lambda - lambda to use in analysis (can be a number or 'auto')
        %   filedir - directory to ouptut saved files (with trailing
        %            slash!)
            n_files = length(file_list);
            for i=1:n_files
                try
                    fprintf('=================================================\n[%d / %d]',i,n_files);
                    fprintf(' LOADING %s\n-------------------------------------------------\n',file_list{i});
                    pf = pffind(file_list{i});
                    STRF.save_pstrf_delta(pf, lambda, filedir, 'bin_width', 10, 'n_lags', 24);
                catch err
                     warning(getReport(err));
                end
            end
        end
        
        function [] = save_pstrf_delta(pf, lambda, filedir, varargin)
        % SAVE_PSTRF_DELTA saves results of pstrf_delta for use in batch analysis
        %
        % INPUTS
        %   pf - p2m file
        %   lambda - number or string 'auto'
        %   filedir - directory to output files to (with trailing slash!)

            pf = PFUtil.removeBadTrials(pf);

            % process lambda param, generate string to be used in file names
            auto_lambda = strcmp(lambda,'auto');
            if strcmp(lambda,'auto')==1
                lambda_str = 'auto';
            else
                lambda_str = sprintf('%.2f',lambda);
            end

            % Run analysis
            S = STRF(pf,'feature','delta','auto_lambda',auto_lambda,varargin{:});

            if auto_lambda
                lambda = S.best_lambda;
            end

            % Save Numerical Results to File
            R = {};
            [R.theta, R.theta_err, R.theta_diff] = S.plotThetaDelta(lambda);
            R.src = S.pf.src;
            R.best_lambda = S.best_lambda;
            R.lambdas = S.lambdas;
            R.lambdas_scores = S.lambdas_scores;
            R.lambdas_scores_err = S.lambdas_scores_err;
            R.lambdas_scores_raw = S.lambdas_scores_raw;
            R.lags = 0:S.BIN_WIDTH:S.BIN_WIDTH*S.N_LAGS;
            [~,exper_name,exper_no] = fileparts(S.pf.src);
            pstrf_result = R;
            save([filedir exper_name exper_no '_pstrf_delta_result_lambda_' lambda_str '.mat'],'pstrf_result');

            % Save Figure to File
            filename = [filedir exper_name exper_no '_pstrf_delta_lambda_' lambda_str '.png'];
            
            % Tries to load the 'printable' export setting
            try
                sdf('printable');
            catch err
                warning('Export preset "printable" undefined, using default export settings');
            end
            
            orient portrait;
            print('-dpng',filename);
            close all;
        end
    end
end