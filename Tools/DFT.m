%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION DFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: A.-L. Robisch
%
% Purpose:  Fourier-transformation of a centered 2D array
%
% Input:   field: 2D array to be transformed
%
% Output:  DFTfield: transformed field
%
%
% Date: 22.01.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function DFTfield = DFT(field)

px=size(field,2);            % number of pixel in x
py=size(field,1);            % number of pixel in y
 
DFTfield=1/sqrt(px*py)*fftshift(fft2(ifftshift(field)));     



