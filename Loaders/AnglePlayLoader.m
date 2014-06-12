classdef AnglePlayLoader < handle
    % ANGLEPLAYLOADER complex cell model
    
    properties
        
        file_list
        image_dir
        oris

    end
    
	properties (Constant)
    	GRIDCURV_PATH = '/lab/stimuli/gridcurv/7';
        ANGLEPLAY_PATH = '/lab/stimuli/curvplay_jackson2/angles_medium';
        MAX_WIDTH = 400;
    end
    
    methods
        function APL = AnglePlayLoader(varargin)
            
            p = inputParser;
            addParameter(p,'imagedir', 0);
            parse(p,varargin{:});
            
            if p.Results.imagedir ~= 0
                APL.image_dir = p.Results.imagedir;
            else
                APL.image_dir = AnglePlayLoader.ANGLEPLAY_PATH;
            end
            
            % orientations to rotate images
            % @todo load from p2m file
            APL.oris = 0:45:315;
            
            image_dir_struct = dir([APL.image_dir '/*.png']);
            APL.file_list = cell(length(image_dir_struct),2);
            for i=1:length(image_dir_struct)
                APL.file_list{i, 1} = image_dir_struct(i).name;
               % GL.file_list{i, 2} = sscanf(GL.file_list{i, 1},'%u-');
            end
            clear image_dir_struct;
        end
        
        function [img] = randomStimulus(APL)
        % RANDOMSTIMULUS loads a random image into memory
        %   RETURN img                  (uint8 matrix) grayscale image
        %   RETURN vector               (uint8 vector) orientations
        
            image_i = randi([1 length(APL.file_list)]);
            vector = APL.file_list{image_i,2};
            path = [APL.image_dir '/' APL.file_list{image_i,1}];
            img = APL.loadImage(path);
        end
        
        function [img] = getByImageNumber(APL, pf, inum, downsampled_size)
        % GETBYIMAGENUMBER loads an image based on the INUM stored
        % in the p2m file, PF.
        % NB: the inum index should be zero-based, since this is how they
        % are stored in the p2m file ev_e entry
        
            % load image info from p2m file,
            image_info = pf.rec(1).params.IMAGE_INFO{inum + 1};
            
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
        
        function img = loadImage(APL, path, rotation, downsampled_size)
            
            try
                [img, map, a] = imread(path);
                % convert indexed png's to grayscale
                % most stimulus png's are not indexed but there are a few...
                if ~isempty(map)
                    img = ind2gray(img,map);
                end
            catch err
                fprintf('path:\n%s\n',path)
                rethrow(err);
            end
            
            img_size = size(img);
            if(length(img_size)==3)
                img = double(rgb2gray(img));
            else
                img = double(img);
            end

            % rescale from [0,255] -> [0,1]
            img = img / 255;
            
            % resize to be MAX_WIDTH
            if size(img,1) > APL.MAX_WIDTH
                img = imresize(img,[APL.MAX_WIDTH, APL.MAX_WIDTH]);
                img(img < 0) = 0; % get rid of scaling artifacts
                img(img > 1) = 1;
            end
            
            % eliminate border effects (some pixels are slightly off-white)
            img(img > 0.9) = 1;
            stim_w = size(img,1);
            
            % invert image since imrotate adds zeros to in blank spaces
            % which in an un-inverted image appear as black. we'll convert
            % it back later...
            img_inverted = imcomplement(img);
            
            % rotate the image and clean up the resulting mess
            % - undo the inversion
            % - the rotation increases the size of the image, and so it
            %   must be cropped back to normal size
            if nargin < 3
                rotation = APL.oris(randi(length(APL.oris),1));
            end
            img_rotated_inverted = imrotate(img_inverted, rotation);
            
            img_rotated = imcomplement(img_rotated_inverted);
            rotated_dim = size(img_rotated,1);
            min_pad = round((rotated_dim - stim_w)/2);
            img_rotated = img_rotated(min_pad+1:(min_pad + stim_w - 1),min_pad+1:(min_pad + stim_w - 1));
            
            % set background to gray, and then demean so that gray is 0
            % effectively rescales range from [0,1] -> [0,0.5] -> [-0.5,0.5]
            img = (img_rotated / 2);% - 0.5;
            
            if nargin >= 4
                slice_w = round(size(img,2) / 2);
                padding = round((size(img,2) - slice_w)/2);
                img = img(padding+1:(padding+slice_w - 1),padding+1:(padding+slice_w - 1));
                img = imresize(img,[downsampled_size, downsampled_size]);
            end
            % @todo change color of curve based on polarity
        end

    end
end