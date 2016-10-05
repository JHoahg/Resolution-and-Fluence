%returns the middle quadrant of size m×n of an image M×N image
%Author J. Hagemann  mJD 57121
%example:
% A = rand(2048);
% figure(); imagesc(A);
% % shows the middle quadrant of matrix A
% B = mid(A, 512, 512);
% figure(); imagesc(B)


function out = mid(in, m, n)
M = size(in, 1);
N = size(in, 2);

if(nargin == 1)
m = floor(M /2) +1;
n = floor(N /2) +1;
end

% JH: i have in my simulations a parameter object which holds some size
% information, thus this is a special case for me, when i just want to pass
% that object
if(nargin == 2)
x = m;
M = size(in, 1);
N = size(in, 2);

m = x.height;
n = x.width;
end


out = in(floor((M-m)/2 + 1) : floor((M-m)/2 + m),...
    floor((N-n)/2 + 1):floor((N-n)/2 + n));
% (p.width2 - p.width)/2 + 1 : (p.width2 - p.width) / 2 + p.width, ...
%       (p.height2 - p.height)/2 + 1 : (p.height2 - p.height) / 2 + p.height));

end