function [result] = RAAR_ff(M, psi, iterations, p)
% inputs:
% M             - measurement (amplitudes)
% psi           - initial guess
% iterations    - number of iterations
% p             - parameters

h = waitbar(0, 'progress');

% get input arrays and put them on GPU
M = gpuArray(M);
%is used in support constraint to ensure energy conservation
% I_tot = sum(M(:));
supp = gpuArray(logical(p.supp));

% Amp_valid can be used to discriminate bad pixels
Amp_valid = gpuArray(p.Amp_valid);
psi = gpuArray(psi);

%  algorithm parameters from data struct p
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
%     sum(abs(psi(:)))
    psi = DFT(psi);
%     sum(abs(psi(:)))
    psi(Amp_valid) = psi(Amp_valid)./ abs(psi(Amp_valid)) .* M(Amp_valid);
%     sum(abs(psi(:)))
%     psi = psi ./ sum(psi(:)) .* I_tot;
    P_M = IDFT(psi);
%     sum(abs(psi(:)))
    
    
    % Reflection on M
    R_M = 2 * P_M - psi_old;
    
%     figure(9000)
%     imagesc(mid(abs(P_M),p)); title('abs pm')
%     figure(9001)
%     imagesc(mid(angle(P_M),p)); title('angle pm')
    
    
    %project S
    % pure phase constraint
    psi(supp) = 1 .* exp(1i .* angle(R_M(supp)));
    % support and negativity constraint
    psi(~supp | angle(psi) > 0) = 1;
%     psi(~supp) = 0;
%     psi( angle(psi) > 0) = 1;
%     psi = psi ./ sum(psi(:)) * I_tot;
    %     mean(abs(psi(:)))
%     figure(9002)
%     imagesc(mid(abs(psi),p)); title('abs ps')
%     figure(9003)
%     imagesc(mid(angle(psi),p)); title('angle ps')
    
    
    % Reflect on S
    psi = 2*psi - R_M;
    
    %  new RAAR iterate
    psi = (b/2) * (psi + psi_old) + (1 -b)*P_M;
%      psi = psi ./ sum(abs(psi(:))) * I_tot;
%      sum(abs(psi(:)))
end
%final projection on M
psi = DFT(psi);
psi(Amp_valid) = psi(Amp_valid)./ abs(psi(Amp_valid)) .* M(Amp_valid);
psi = IDFT(psi);

result = gather(psi);
close(h);
end
