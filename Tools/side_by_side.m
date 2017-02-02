% makes side by side plots of images
% Parameters:
% im1 : image for the left half
% im2: image for the right half
% upper & lower: lower/upper percent of plot values to be included in colormap scaling.
% JH 57605 mJD 
function side_by_side(im1, im2, lower, upper)
    if(nargin < 4)
        lower = 0;
        upper = 100;
    end

    im1 = gather(im1);
    im2 = gather(im2);
    [M,N] = size(im1);
    
    if(size(im1) ~= size(im2))
        error('images have unequal size')
    end
    
    new_img = zeros(M,N);
    
    new_img(:, 1:floor(N/2)) = im1(:, 1:floor(N/2));
    new_img(:, ceil(N/2):end) = im2(:, ceil(N/2):end);
    figure
    imagesc(new_img)
    %adapt contrast to 5-95%
    range = prctile(new_img(:), [lower upper]);
    caxis([range(1) range(2)])
    
    axis_pos = axis;
    posx = axis_pos(2)*0.5;
%     posy = axis_pos(4)*0.5;
    
    line([posx posx],[0 M],'LineWidth',1,'Color','r');


end