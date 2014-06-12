function theta = pstrf(pf, lambda, feature, varargin)
% PSTRF plots the strf for the given PF file

if nargin < 2
    lambda = 0;
end

if nargin < 3
    feature = 'binary';
end

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

S.plotTheta(lambda);

% subplot(3,3,[6,9]);
% if strcmp(feature,'image')
%     ylabel('Stimulus Features (Pixels)');
% end

end

