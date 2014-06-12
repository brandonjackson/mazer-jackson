function [kernel,pos] = gabor_filter_translate(gabor, stimulus_size, x_offset, y_offset)
%GABORFILTERTRANSLATE moves a gabor filter to a position within a kernel the
%size of the stimulus.

gabor_size = size(gabor);

kernel = zeros(stimulus_size,'single');

% If stimulus larger than kernel, add appropriate padding to gabor filter to fit
if stimulus_size(1) >= gabor_size(1)
    origin = stimulus_size/2;
    center = round([(origin(1)+x_offset), (origin(2)-y_offset)]); 
    % note: y_offset subtracted since y axis is flipped in matlab

    left = center(2) - floor(gabor_size(2)/2);
    right = center(2) + floor(gabor_size(2)/2);
    top = center(1) - floor(gabor_size(1)/2);
    bottom = center(1) + floor(gabor_size(1)/2);

    pos = [left right top bottom];

    kernel(left:right,top:bottom) = gabor;

% If stimulus smaller than kernel, crop the gabor filter to fit
else
    overflow = round((gabor_size(1) - stimulus_size(1)) / 2);
    kernel(:,:) = gabor((overflow+1):(overflow + stimulus_size(1)),(overflow+1):(overflow + stimulus_size(1)));
    pos = [1 stimulus_size(1) 1 stimulus_size(2)];
end

end