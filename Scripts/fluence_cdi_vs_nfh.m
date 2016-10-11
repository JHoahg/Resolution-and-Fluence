%% play with dose in CDI / Holography
% change folder to Scripts before executing, sadly one can not find out
% in Matlab where a plain script lies in the filesysytem :(
% http://de.mathworks.com/matlabcentral/newsreader/view_thread/269433
working_dir = '.';
cd(working_dir)
%data dir
TB_path = '../Tools';
addpath(genpath(TB_path));

path = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH'); %clears path
path = [path ':/usr/local/cuda/lib'] % new path
setenv('LD_LIBRARY_PATH', path);
clear path
% 
set(0,'DefaultFigureColor','w')
cmap = flipud(colormap('bone(512)'));
set(0,'DefaultFigureColormap', cmap);

cmd = sprintf(['axis image; colorbar;  set(gcf,''PaperPosition'',[0 0 5 5], ''PaperSize'', [5 5]);',...
    ' set(gca,''xtick'',[],''ytick'',[])']);
set(0,'DefaultImageCreateFcn', cmd)
set(groot,'defaultLineLineWidth',1.5)
clear cmap cmd;
%% Simulation parameters
p.width = 512;
p.height = 512;
p.width2 = 1024;
p.height2 = 1024;
p.rec_width = 1024;
p.rec_height = 1024;

p.F = 1e-3;iaepns = 200;
p.b_0 = 0.99;
p.b_m = 0.75;
p.b_s = 150;
p.num_iterations = 200;
p.Amp_valid = logical(true(p.rec_height, p.rec_width));
%% Sample Setup
if(1) %dicty sketch
    sample = prepare_probe('dicty_sketch.png', 'dicty_sketch.png',  0, 1, 1, 1, p);
    sample = abs(sample) .* exp(1i.*(angle(conj(sample))));
    p.supp = load('support_dicty_ff.mat');
    p.supp = p.supp.supp.I_supp;
    p.supp = imresize(p.supp, (p.height/256), 'box');
end

if(0) % use binary bitmap
    upscale =  10;
    sample = binary_bitmap(10, 10, -1, upscale, p.rec_height, p.rec_width);
    p.supp = ones(10*upscale, 10*upscale);
    p.supp(1:floor(10*upscale/2), 1:ceil(10*upscale/2)) = 0;
end

p.supp = pad_to_size(p.supp, p.rec_height, p.rec_width, 'zero');

figure()
imagesc(angle(mid(sample,p)));

figure()
imagesc( mid(p.supp, p));


%% holography
prop = Propagator(p.F, p.F, p.width2, p.height2,1);
id_holo = abs(prop.propTF(sample)).^2;
id_holo = id_holo ./ sum(id_holo(:));
figure; imagesc(id_holo)
%%
holo = imnoise((id_holo*p.num_photons.* numel(id_holo))*1e-12, 'poisson')*1e12;
holo = holo ./ sum(holo(:)) .* numel(holo);
if(1) % use RAAR
    holo_rec{1} = RAAR_nf(sqrt(double(holo)), ones(p.rec_height), p.num_iterations,...
        p.F, p);
    
end

if(1)
    holo = imnoise((id_holo*p.num_photons.* numel(id_holo))*1e-12, 'poisson')*1e12;
    holo = holo ./ sum(holo(:)) .* numel(holo);
    
    holo_rec{2} = RAAR_nf(sqrt(holo), ones(p.rec_height),...
        p.num_iterations, p.F, p);
end

if(1) %use HIO
    hio_settings = hio_fresnel();
    hio_settings.use_gpu = 1;
    hio_settings.purephase_assumption = 1;
    hio_settings.supp = p.supp;
    hio_settings.N_it = p.num_iterations;
    [holo_rec{2}, ~] = hio_fresnel(holo, p.F, p.F, hio_settings);
end


%% compare with frc
p.beta = 7;
frc_holo = FSC(angle((holo_rec{1})),...
    angle((holo_rec{2})),p );


frc_holo_to_ori = FSC(angle((sample)),...
    angle((holo_rec{1})), p);


%% CDI
id_FT = abs(DFT(sample)).^2;
id_FT = id_FT ./ sum(id_FT(:));
figure; imagesc(log10(abs(id_FT)));
%%
FT = imnoise((id_FT*p.num_photons.* numel(id_FT))*1e-12, 'poisson')*1e12;
FT(FT == 0) = 2*eps;
FT = FT ./ sum(FT(:)) .* numel(FT);
%%

if(1)
    tic
    [CDI_rec{1}] = RAAR_ff(sqrt(FT),...
        exp(1i .* rand(p.rec_height)), p.num_iterations, p);
    toc
end

if(1)
    FT = imnoise((id_FT*p.num_photons.* numel(id_FT))*1e-12, 'poisson')*1e12;
    FT(FT == 0) = 2*eps;
    FT = FT ./ sum(FT(:)) .* numel(FT);
    
    [CDI_rec{2}] = RAAR_ff(sqrt(FT),...
        exp(1i .* rand(p.rec_height)), p.num_iterations, p);
end
if(0) %use HIO
    hio_settings_ff = hio_farfield();
    hio_settings_ff.use_gpu = 1;
    hio_settings_ff.purephase_assumption = 1;
    hio_settings_ff.supp = p.supp;
    hio_settings_ff.N_it = p.num_iterations;
    [CDI_rec{2}, result] = hio_farfield(FT, hio_settings_ff);
end

%% compare with frc
p.beta = 7;
[err, greg] = dftregistration(fft2((CDI_rec{1})), fft2((CDI_rec{2})), 3);
frc_CDI = FSC(angle((CDI_rec{1})), angle(ifft2(greg)),p);


[err, aligned_sample] = dftregistration(fft2(angle(CDI_rec{1})),...
    fft2(angle(sample)), 5);
%             err
aligned_sample = (ifft2(aligned_sample));
frc_CDI_to_ori = FSC(angle((CDI_rec{1})), aligned_sample, p);

%% compare both
% two recons
figure
plot(frc_CDI.nu, abs(frc_CDI.frc), 'Displayname', 'CDI')
hold on
plot(frc_holo.nu, abs(frc_holo.frc), 'DisplayName', 'holo')
plot(frc_CDI.nu, abs(frc_CDI.T_hbit))
hold off
legend('show')
title(sprintf('FRC CDI %g photons', p.num_photons))
xlim([-0.01 0.5])
%% to ori.
figure
plot(frc_CDI_to_ori.nu, abs(frc_CDI_to_ori.frc), 'Displayname', 'CDI')
hold on
plot(frc_holo_to_ori.nu, abs(frc_holo_to_ori.frc), 'DisplayName', 'holo')
plot(frc_CDI_to_ori.nu, abs(frc_CDI_to_ori.T_hbit))
hold off
legend('show')
title(sprintf('FRC CDI to original %g photons', p.num_photons))
xlim([-0.01 0.5])
ylim([0 1.05])
%% images...
figure;
imagesc(angle(mid(holo_rec{1},p)))
title('Result of holographic reconstruction')
c = colorbar; c.Label.String = 'Phase (rad)';

figure;
imagesc(angle(mid(CDI_rec{1},p)))
title('Result of CDI reconstruction')
c = colorbar; c.Label.String = 'Phase (rad)';

