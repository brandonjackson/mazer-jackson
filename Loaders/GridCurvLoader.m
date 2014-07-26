classdef GridCurvLoader < SuperLoader
    % GRIDCURVLOADER loads gridcurv stimuli
    
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
        %   RETURN img                  (uint8 matrix) grayscale image
        %   RETURN vector               (uint8 vector) orientations
        
            image_i = randi([1 length(GCL.file_list)]);
            vector = GCL.file_list{image_i,2};
            path = [GCL.image_dir '/' GCL.file_list{image_i,1}];
            img = GCL.loadImage(path);
        end
        
        function img = loadImage(GCL, path)
            try
                img = imread(path);
            catch err
                fprintf('path:\n%s\n',path)
                rethrow(err);
            end

            img = uint8(rgb2gray(img));
        end

    end
end