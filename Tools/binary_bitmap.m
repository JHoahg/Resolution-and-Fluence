% creates binary bitmap image á la Jahn et al
% Author: JH,  mJD 57595
%
% parameters:
% nx & ny : dimensions of random binary bitmap image
% shift : phase shift for pixels == 1
% upscale : blow up binary pixels e.g. 1x1 -> 3x3
% width, height : pad binary image to given width and height

% example:
% bitmap = binary_bitmap(10, 10, -0.25, 3, 512, 512);
%figure; imagesc(angle(bitmap))
function [sample] = binary_bitmap(ny, nx, shift, upscale, width, height)
sample = rand(nx * ny, 1);
sample(sample >= 0.5) = 1;
sample(sample < 0.5) = 0;

sample = reshape(sample, ny, nx);
% break up symmetry
sample(1:nx/2, 1:ny/2) = 0;

sample = imresize(sample, upscale, 'box');
sample = pad_to_size(sample, width, height, 'zero');
sample = exp(1i * sample * shift);
end