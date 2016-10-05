% prepare probes
% this function takes 2 paths as arguments which are intaerpreted as
% amplitude and phase of a wavefield(probe or object or whatever). the user
% can set the ranges of amplitude and phase. p is a parameter struture, its
% purpose here is to give the size of the desired wavefield.
% JH 20141204
function [beam] = prepare_probe(phapath, amppath, lower_phase, upper_phase,...
                    lower_amp, upper_amp, p, gauss_filt_fwhm)

if nargin < 8
    gauss_filt_fwhm = 0;
end
                
% read in and convert to grayscale
amp = imread(amppath);

% convert only to grayscale if we hav a color image, otherwise rgb2gray
% throws an error
if ndims(amp) > 2
   amp = rgb2gray(amp);
end
   
pha = imread(phapath);

if ndims(pha) > 2
    pha = rgb2gray(pha);
end

amp = imresize(amp, [p.width p.width]);
pha = imresize(pha, [p.width p.width]);

if gauss_filt_fwhm > 1e-6
    amp = imgaussfilt(amp, gauss_filt_fwhm/(2*sqrt(2*log(2))));
    pha = imgaussfilt(pha , gauss_filt_fwhm/(2*sqrt(2*log(2))));
end

% scaled amplitudes
amp_range = double(max(amp(:)) - min(amp(:)));
amp = lower_amp + (double(amp) ./ (amp_range)) .* (upper_amp - lower_amp);

% scaled phases
pha_range = double(max(pha(:)) - min(pha(:)));
pha = lower_phase + (double(pha) ./ (pha_range)) .* (upper_phase - lower_phase);

beam = amp .* exp(1i .* pha);
beam = pad_fadeout(beam, [p.height2, p.width2],(p.width2-p.width)*0.01, exp(1i*0));
end