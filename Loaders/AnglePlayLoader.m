classdef AnglePlayLoader < SuperLoader
    % ANGLEPLAYLOADER load AnglePlay stimuli
    %
    % @todo implement getByTrigger()
    % @todo implement getByStimulusParams()
    
    properties
        
        % Directory where images are stored
        image_dir
        
        % List of image files
        image_list
        
        % List of orientations that images are rotated by
        oris

        % Task Version (0==old, 1==new)
        task_version
    end
    
	properties (Constant)
        MAX_WIDTH = 400;
    end
    
    methods
        function APL = AnglePlayLoader(pf, varargin)
            
            APL = APL@SuperLoader(pf); % Construct superclass
            
            APL.image_dir = APL.pf.rec(1).params.imagedir;
            
            % Load list of orientations to rotate images
            oristep = APL.pf.rec(1).params.oristep;
            APL.oris = 0:oristep:(360-oristep);
            
            % Load list of file names
            image_dir_struct = dir([APL.image_dir '/*.png']);
            APL.image_list = cell(length(image_dir_struct),1);
            for i=1:length(image_dir_struct)
                APL.image_list{i, 1} = image_dir_struct(i).name;
            end
            clear image_dir_struct;
            
            % Determine whether old or new version of task used
            APL.task_version = AnglePlayUtil.taskVersion(APL.pf);
        end
        
        function [img] = randomStimulus(APL)
        % RANDOMSTIMULUS loads a random image into memory
        %   RETURN img                  (uint8 matrix) grayscale image
        %   RETURN vector               (uint8 vector) orientations
        
            image_i = randi([1 length(APL.image_list)]);
            path = [APL.image_dir '/' APL.image_list{image_i,1}];
            rotation = APL.oris(randi(length(APL.oris),1));
            img = APL.loadImage(path, rotation);
        end
        
        function [img] = getByImageNumber(APL, inum, downsampled_size)
        % GETBYIMAGENUMBER loads an image based on the INUM stored
        % in the p2m file, PF.
        % NB: the inum index should be zero-based, since this is how they
        % are stored in the p2m file ev_e entry
        
            % load image info from p2m file,
            image_info = APL.pf.rec(1).params.IMAGE_INFO{inum + 1};
            
            % tokenize image info, extracting rotation and filename
            tokens = strsplit(image_info,[char(9) ' ']); % split by tab+space
            filename = tokens{1};
            path = [APL.image_dir '/' filename];
            rotation = str2num(tokens{5});
            
            if nargin < 4
                img = APL.loadImage(path,rotation);
            else
                img = APL.loadImage(path,rotation,downsampled_size);
            end
        end
        
%         function img = getByTrigger(APL, pf, trigger)
%             
%         end
        
        function stimulus = loadImage(APL, path, rotation)
            
            try
                [img, map, a] = imread(path);
                % convert indexed png's to grayscale
                % most stimulus png's are not indexed but there are a few...
                if ~isempty(map)
                    img = ind2gray(img,map);
                end
                
                % Convert RGB to grayscale
                if(length(size(img))==3)
                    img = rgb2gray(img);
                end
                
            catch err
                fprintf('path:\n%s\n',path)
                rethrow(err);
            end
            
            % Convert to double, rescale from [0,255] -> [0,1]
            img = im2double(img); 
            
            % Eliminate border effects (some pixels are slightly off-white)
            img(img > 0.9) = 1;
            
            % Do Some Math First!
            % Calculate image size (size to scale the image to)
            if APL.task_version==1
                image_size = (10 * APL.pf.rec(1).params.rfsigma) * APL.pf.rec(1).params.ds;
            elseif isfield(APL.pf.rec(1).params,'scale')
                image_size = size(img,1) * APL.pf.rec(1).params.scale;
            else
                image_size = size(img,1);
            end
            
            % Calculate stimulus size (size of final stimulus (pre-scaling))
            stimulus_size = round(image_size * 1.5); % in theory, should be * sqrt(2) 
            
            % Resize the Image
            img = imresize(img, [image_size image_size]);
            img(img < 0) = 0; % get rid of scaling artifacts
            img(img > 1) = 1;
            
            % Rotate the Image
            % invert image (since imrotate adds zeros to in blank spaces
            % which in an un-inverted image appear as black), rotate, and
            % un-invert back to normal color range
            img_inverted = imcomplement(img);
            img_rotated_inverted = imrotate(img_inverted, rotation);
            img_rotated = imcomplement(img_rotated_inverted);
            
            % Insert rotated image into big stimulus foreground frame
            % (with padding to account for changes in size when rotated)
            rotated_dim = size(img_rotated,1);
            pad = round((stimulus_size - rotated_dim) / 2);
            stimulus_fg = ones([stimulus_size,stimulus_size]);
            stimulus_fg(pad:(pad + rotated_dim - 1),pad:(pad + rotated_dim - 1)) = img_rotated;
                        
            % Set background to gray, rescale from [0,1] -> [0,0.5]
            stimulus_fg = (stimulus_fg / 2);
            
            % Create Gaussian Envelope, rescale so range is [0,1]
            if APL.task_version==0
                sigma = APL.pf.rec(1).params.sigma;
            else
                sigma = APL.pf.rec(1).params.rfsigma * APL.pf.rec(1).params.nsigma;
            end
            gaussian_envelope = fspecial('Gaussian',stimulus_size,sigma);
            gaussian_envelope = gaussian_envelope / max(gaussian_envelope(:));
            
            % Create Gray Stimulus Background
            stimulus_bg = 0.5 * ones(size(stimulus_fg));
            
            % Blend Stimulus Foreground + Alpha with Background
            stimulus = gaussian_envelope.*stimulus_fg + (1-gaussian_envelope).*stimulus_bg;
            
            % @todo change color of curve based on polarity
        end

    end
end