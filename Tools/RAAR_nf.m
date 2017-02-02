function [result] = RAAR_nf(M, psi, iterations, F, p)
% inputs:
% M             - measurement (amplitudes)
% psi           - initial guess
% iterations    - number of iterations
% F             - Fresnel number for propagation
% p             - parameters

h = waitbar(0, 'progress');

% get input arrays and put them on GPU
M = gpuArray(M);
supp = gpuArray(logical(p.supp));

% Amp_valid can be used to discriminate bad pixels
Amp_valid = gpuArray(p.Amp_valid);
psi = gpuArray(psi);

% Propagator objects
if(isfield(p,'oversample') == 0)
    warning('using default oversampling factor = 1.');
    oversample = 1;
else
    oversample = p.oversample;
end
prop = PropagatorGPU(F, F, size(M,2), size(M,1), oversample);
prop_back = PropagatorGPU(-F, -F, size(M,2), size(M,1), oversample);

b_0 = p.b_0;
b_m = p.b_m;
b_s = p.b_s;

for ii = 1:iterations
    waitbar(ii / iterations, h, ...
        sprintf('%d / %d',ii, iterations));


% relaxation parameter for current iteration
    b = exp(-(ii/b_s)^3)*b_0 + (1 - exp(-(ii/b_s)^3))*b_m;
    
    psi_old = psi;
    
    % P_M
    psi = prop.propTF(psi);
    psi(Amp_valid) = psi(Amp_valid)./ abs(psi(Amp_valid)) .* M(Amp_valid);
    psi = psi./ abs(psi) .* M;
    P_M = prop_back.propTF(psi);
    
    % Reflection on M
    R_M = 2 * P_M - psi_old;
    
    %project S
    % pure phase constraint
    psi(supp) = exp(1i .* angle(R_M(supp)));
    % support and negativity constraint
    psi(~supp | angle(psi) < 0) = 1;
        
    % Reflect on S
    psi = 2*psi - R_M;
    
    %  new RAAR iterate
    psi = (b/2) * (psi + psi_old) + (1 -b)*P_M;

end
% P_M
    psi = prop.propTF(psi);
    psi(Amp_valid) = psi(Amp_valid)./ abs(psi(Amp_valid)) .* M(Amp_valid);psi = psi./ abs(psi) .* M;
    psi = prop_back.propTF(psi);

    result = gather(psi);
    close(h);
end