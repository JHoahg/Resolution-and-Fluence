%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION IDFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: A.-L. Robisch
%
% Purpose:  Fourier-back-transformation of a centered 2D array
%
% Input:   field: 2D array to be transformed
%
% Output:  IDFTfield: transformed field
%
%
% Date: 22.01.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function IDFTfield = IDFT(field)

        px=size(field,2);
        py=size(field,1);
        IDFTfield=sqrt(px*py)*fftshift(ifft2(ifftshift(field)));

