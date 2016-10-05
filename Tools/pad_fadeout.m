%=========================================================================
% Pads a 2D input image (or a stack of such) by a smooth transition from
% replicate-padding to a constant value, i.e. a fade out of the image's
% boundary values
%
% Example:
%       im = rand([64,64,3]);
%       res = pad_fadeout(im, [128,128], 10, 1);
%       imagesc(res(:,:,2));
%
% Last modified on June, 9th 2015 by Simon Maretzke 
%=========================================================================


function res = pad_fadeout(im_in, N_out, transition_length, padval, method)

% Input arguments:
%
%            im_in1: N_in x X array containing the size N_in input images
%             N_out: 2-tuple containing the desired size of the padded output images.
% transition_length: Length of the transition region (in pixels) between
%                    replicate- and constant padding (optional, default = 100).
%            padval: Constant padding value to be applied beyond the replicate
%                    region (optional, default = 0)
%
% Output arguments:
%       
%               res: N_out x X array containing the padded output images

        % Assign default values
        if nargin < 5
                method = 'replicate';
        end
        if nargin < 4
                padval = 1.;
        end
        if nargin < 3
                transition_length = 100;
        end

        N_in = size(im_in);
        pad_low = ceil( (N_out - N_in(1:2)) / 2 );      % Number of lower and
        pad_up = floor( (N_out - N_in(1:2)) / 2 );      % upper pixels to be padded

        
        % Construct mask that is 1 within the region covered by the original
        % image and decays like a cosine to 0 in the transition region.
        [X,Y] = ndgrid( (pi/transition_length) * min( transition_length, max( 0, abs( (0:N_out(1)-1).' - (pad_low(1)+N_in(1)/2) ) - N_in(1)/2 ) ), ...
                        (pi/transition_length) * min( transition_length, max( 0, abs( (0:N_out(2)-1).' - (pad_low(2)+N_in(2)/2) ) - N_in(2)/2 ) )      );        
        transition_mask = 0.25 * ( 1 + cos(X) ) .* ( 1 + cos(Y) );

        % Superimpose constant and replicate-padded image with the transition mask.
        res = padval * ones([N_out, N_in(3:end)]);
        res = res + bsxfun(@times, transition_mask, padarray( padarray(im_in, pad_low, method, 'pre'), pad_up, method, 'post') - res );

end
