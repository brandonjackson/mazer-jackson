function grk2(file, lat, winsize, varargin)
%function grk(file, lat, winsize, varargin)
%
% plot full gratrev kernel (grk) at fixed lag
%
%    grk(file, lat, winsize, [options])
%
%  INPUT
%    file      - pf struct, p2m filename etc..
%    lat       - temporal latency (ms)
%    winsize   - temporal integration window (ms)
%
%    options - specified as <'option', value> pairs:
%      smooth    - 'spatial' smoothing for 2D kernel plot
%      tsigma    - temporal smoothing width in ms psth
%      cpd       - 0/1 plot in cycles/deg instead of cycles/pix -- p2m only!
%      unit      - unit descriptor (see p2mgetspikes)
%                  nb: uniselect has precedence
%

if nargin < 3
  error('specify at least pf/file, lat and winsize -- see help');
end

dolog = 1;

opt.lat = lat;
opt.winsize = winsize;

opt.smooth = 0;                         % 'spatial' smoothing
opt.tsigma = 5;				% PSTH temporal smoothing window (ms)
opt.unit = 'ttl';			% default to ttl
opt.cpd = 0;                            % use cycles/deg?
opt = getopts(opt, varargin{:});

pf = getpf(file);
% See if it's a non-cart gratrev file..
if ~isempty(strfind(pf.rec(1).params.X_version_info{2}, 'ncgratrev'))
  nc = 1;
elseif ~isempty(strfind(pf.rec(1).params.X_version_info{2}, 'gratrev'))
  nc = 0;
else
  error('not gratrev file');
end
  
if isfield(pf, 'uniselect')
  code = ['UNI:' pf.uniselect];
elseif ~isempty(opt.unit)
  % pre-unispike
  p2mgetspikes(opt.unit);
  code = p2mgetspikes;
else
  code = p2mgetspikes;
end

file = pf.src;
code = p2mgetspikes;

respmat = [];
for n = 1:length(pf.rec)
  [ix, ts] = p2mFindEvents(pf, n, 'FLIP');
  for k=1:length(ix)
    s = pf.rec(n).ev_e{ix(k)};
    if ~nc
      s = str2num(s(1+length('FLIP'):end));
      if length(s) < 4
        s(4) = 1.0;
      end
      ori = s(1); sf = s(2); pha = s(3); cont = s(4);
      sf = sf / (2 * pf.rec(1).params.radius + pf.rec(1).params.taper);
    else
      s = sscanf(s, 'FLIP s %f %f');
      if isempty(s)
        % not a sine-grating, skip it..
        continue;
      end
      ori = s(1); sf = s(2); pha = 0; cont = 1.0;
      sf = sf / (2 * pf.rec(1).params.radius + pf.rec(1).params.taper);
    end
    if opt.cpd
      sf = round(100 * sf * pf.rec(n).params.mon_ppd) / 100.0;
    end
    
    % accumulate response matrix
    nspikes = length(p2mgetspikes(pf, n, ...
                                  ts(k)+opt.lat, ts(k)+opt.lat+opt.winsize));
    spikes_per_sec = nspikes / opt.winsize * 1000.0;
    respmat = [respmat; spikes_per_sec ori sf pha cont];
  end
end

pats = {{2, 3, 'ori'} {3, 2, 'sf'}};
alist = [];
for k = 1:length(pats)
 
  %subplot(1+length(pats), 2, [k*2-1 k*2]);
  alist = [alist gca];
  
  hlist = [];
  tlist = {};
  
  a = pats{k}{1};
  b = pats{k}{2};
  
  v1s = unique(respmat(:,a))';
  v2s = unique(respmat(:,b))';
  cm = jet(1+max(length(v1s),length(v2s)));
  n = 1;
  m = [];
  for s = v2s
    x = [];  y = []; e = [];
    for o = v1s
      ix = find(respmat(:,a) == o & respmat(:,b) == s);
      x = [x o];
      y = [y mean(respmat(ix, 1))];
      e = [e sem(respmat(ix, 1))];
    end
    m = [m; y];
   
    n = n + 1;
  end
end

%subplot(1+length(pats), 2, (1+length(pats))*2-1)
if opt.smooth
m = smooth2d(v2s, v1s, m, opt.smooth);
end
%imagesc(v2s, v1s, m');
%hold on;
[C,h]=contourf(v2s, v1s, m',50);

try
    shading flat;
catch err
end

set(gca, 'ydir', 'normal');

if dolog, set(gca, 'YScale', 'log'); end
xlabel(pats{1}{3});
ylabel(pats{2}{3});
axis tight;
colorbar;

function m = smooth2d(xs, ys, m, smooth)

[X, Y] = meshgrid(xs, ys);
X = X - mean(X(:));
Y = Y - mean(Y(:));
xsigma = smooth;
ysigma = smooth;
sk = exp(-(((X).^2))/(2*xsigma.^2)) .* exp(-(((Y).^2))/(2*ysigma.^2));
sk = sk ./ sum(sk(:));
if ~all(all(isnan(sk)))
  m = conv2(m, sk, 'same');
end
  
  
