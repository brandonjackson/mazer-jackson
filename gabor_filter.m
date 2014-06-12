function gb = gabor_filter(params, varargin)
%GABORFILTER Create a Gabor filter.
%   GB = GABORFILTER(THETA, WAVELENGTH, PHASE, SIGMA, ASPECT, KERNEL_SIZE)
%   creates a Gabor filter GB with orientation THETA (in radians),
%   wavelength WAVELENGTH (in pixels), phase offset PHASE (in radians),
%   envelope standard deviation SIGMA, aspect ratio ASPECT, and dimensions
%   KERNEL_SIZE x KERNEL_SIZE. KERNEL_SIZE is an optional parameter, and if omitted default
%   dimensions are selected.

p = inputParser;
addRequired(p,'params');
addParamValue(p,'phase',10000,@isnumeric);
addParamValue(p,'theta',10000,@isnumeric);
parse(p,params,varargin{:});

% override default phase
if p.Results.phase ~= 10000
    params.phase = p.Results.phase;
else
    params.phase = 0;
end

% override default theta
if p.Results.theta ~= 10000
    params.theta = p.Results.theta;
end

theta = degtorad(params.theta + 90); % convert to radians!
wavelength = params.wavelength;
phase = params.phase;
sigma = params.sigma;
aspect = params.aspect;
if ~isfield(params,'kernel_size')
    kernel_size = 8*sigma*aspect;
else
    kernel_size = params.kernel_size;
end

if numel(kernel_size) == 1
  kernel_size = [kernel_size kernel_size];
end

xmax = floor(kernel_size(2)/2);
xmin = -xmax;
ymax = floor(kernel_size(1)/2);
ymin = -ymax;

[xs, ys] = meshgrid(xmin:xmax, ymax:-1:ymin);

%% My Code

% Empty gabor filter matrix
gb = zeros((xmax-xmin),(ymax-ymin));

% Translate x and y values into x' and y'
xs_theta = cos(theta)*xs + sin(theta)*ys;
ys_theta = -sin(theta)*xs + cos(theta)*ys;

% Gabor Matrix
gb = sin((2*pi/wavelength).*ys_theta + phase).*exp(-1*((xs_theta.^2/aspect^2)+ys_theta.^2)/(2*sigma^2));

% Normalize to 0 by subtracting mean
gb = gb-mean(gb(:));

% Scale max response to 1
% gb = gb / max(gb(:));
% max(gb(:))

% Convert to Single
gb = single(gb);