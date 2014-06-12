function theta = pstrf_delta(pf, lambda,  varargin)
% PSTRF plots the strf for the given PF file using derivative
% stimulus representation (i.e. stims as delta functions)
% i.e. instead of a 5 ms stimulus represented as 0 1 1 1 1 1 0
%                           it is represented as 0 1 0 0 0 0 0

if nargin < 2
    lambda = 0;
end

feature = 'delta';

% if nargin < 4
%     jackknife = 0;
% end

pf = PFUtil.removeBadTrials(pf);

auto_lambda = strcmp(lambda,'auto');

S = STRF(pf,'feature',feature,'auto_lambda',auto_lambda,varargin{:});

% if jackknife > 0
%     theta = PFUtil.jackknife(pf, @STRF.fullRegression, jackknife, 'lambda',lambda, 'feature',feature);
%     jackknife_str = [', ' num2str(jackknife) '-draw jackknife '];
% else
%     theta = S.STRF.fullRegression(pf,'lambda',lambda, 'feature',feature);
%     jackknife_str = '';
% end
jackknife_str = '';

if auto_lambda
    lambda = S.best_lambda;
end

S.plotThetaDelta(lambda, [' STRF (\lambda = ' num2str(lambda) ', feature=' feature jackknife_str ' )']);

% subplot(3,3,[6,9]);
% if strcmp(feature,'image')
%     ylabel('Stimulus Features (Pixels)');
% end

end

