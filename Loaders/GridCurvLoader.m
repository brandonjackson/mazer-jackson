classdef GridCurvLoader < SuperLoader
    % GRIDCURVLOADER loads gridcurv stimuli
    % 
    % WARNING: This class does not actually generate stimuli... it loads
    % random images from a directory that contains pre-generated images.
    % Since the set of possible stimuli is practically infinite, it is
    % impossible to load a specific image from the set of pre-compiled
    % images since it most likely does not exist
    %
    % @todo Rewrite to generate stimulus images
    
    properties
        file_list
        image_dir
    end
    
	properties (Constant)
    	GRIDCURV_PATH = '/lab/stimuli/gridcurv/7';
    end
    
    methods
        function GCL = GridCurvLoader(pf)
            
            GCL = GCL@SuperLoader(pf); % Construct superclass
            
%             if ~isempty(varargin.image_dir)
%                 GCL.image_dir = varargin.image_dir;
%             end
            GCL.image_dir = GridCurvLoader.GRIDCURV_PATH;
            
            image_dir_struct = dir([GCL.image_dir '/*.png']);
            GCL.file_list = cell(length(image_dir_struct),2);
            for i=1:length(image_dir_struct)
                GCL.file_list{i, 1} = image_dir_struct(i).name;
                GCL.file_list{i, 2} = sscanf(GCL.file_list{i, 1},'%u-');
            end
            clear image_dir_struct;
        end
        
        function [img, vector] = randomStimulus(GCL)
        % RANDOMSTIMULUS loads a random wavelet into memory
        %   RETURN img                  (double matrix) image, range [0,1]
        %   RETURN vector               (uint8 vector) orientations
        
            image_i = randi([1 length(GCL.file_list)]);
            vector = GCL.file_list{image_i,2};
            path = [GCL.image_dir '/' GCL.file_list{image_i,1}];
            img = GCL.loadImage(path);
        end
        
        function img = loadImage(GCL, path)
        % LOADIMAGE loads the image of a grid of wavelets
        %   RETURN img                  (double matrix) image, range [0,1]

            try
                img = imread(path);
            catch err
                fprintf('path:\n%s\n',path)
                rethrow(err);
            end

            % Convert image to range [0,1]
            img = im2double(img);
        end

    end
end