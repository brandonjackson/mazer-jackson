classdef GaborUtil < handle
    % GABORUTIL library of static methods for creating/editing gabor
    % filters
    
    properties
    end
    
	properties (Constant)
    end
    
    methods
    end
    
    methods (Static)
        
        function gb = createFilter(ori, wavelength, sigma, phase, aspect)
            
        % CREATEFILTER Create a Gabor filter.
        %   creates a Gabor filter GB with orientation ORI (in radians),
        %   wavelength SF (in pixels), phase offset PHASE (in radians),
        %   envelope standard deviation SIGMA, aspect ratio ASPECT, and dimensions
        %   KERNEL_SIZE x KERNEL_SIZE. KERNEL_SIZE is an optional parameter, and if omitted default
        %   dimensions are selected.
        
            if nargin < 5
                aspect = 1;
            end
            
            if nargin < 4
                phase = 0;
            end

            theta = degtorad(ori + 90); % convert to radians!
            kernel_size = [8*sigma*aspect 8*sigma*aspect];

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

        end
        
    end
end