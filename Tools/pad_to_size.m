% pad_to_size padds array via replication to given size.
% mode = zero padds with zeros
% other modes are 'circular' 'symmetric', reproduce std. matlab behavior
% Author JH 2015012
% example:
% bla = rand(11);
% bla2 = pad_to_size(bla, 25, 25);

function IM = pad_to_size(im, ny, nx, mode)
if nargin ~= 4;
    mode = 'replicate';
end

im_size = size(im);

padx = floor((nx - im_size(2))/2);
pady = floor((ny - im_size(1))/2);

% 'rough' padding
if ~strcmp(mode, 'zero')
    IM = padarray(im, [pady padx], mode);
else
    IM = padarray(im, [pady padx]);
end

% take care of single pixel aberrations
if logical(mod(im_size(2), 2)) || logical(mod(im_size(1), 2))
    IM_size = size(IM);
    
    if ~strcmp(mode, 'zero')
        IM = padarray(IM, [ny-IM_size(1), nx-IM_size(2)],'pre', mode);
    else
        IM = padarray(IM, [ny-IM_size(1), nx-IM_size(2)],'pre');
    end
end
end