function [res] = circshift_sp(A, shift)
% Subpixel accuracy circshift using linear interpolation in real space
% aru 2014

    a=floor(shift(2));
    b=floor(shift(1));
    wx = shift(2)-a;
    wy = shift(1)-b;

    res = circshift(A,[b,a]).*(1-wy) + circshift(A,[b+1,a]).*wy;
    res = res.*(1-wx) + circshift(res,[0,1]).*wx;
end
