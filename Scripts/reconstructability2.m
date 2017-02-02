%% play with dose in CDI / Holography
working_dir = '/homegroups/AG_Salditt/Publikationen/in_preparation_submitted_2016/Hagemann_Salditt_Noise_and_Resolution/reconstructability';
cd(working_dir)
%data dir
TB_path = '/home/AG_Salditt/Projects_X-ray_Imaging/Toolbox/release/';

% addpath('/home/jhagemann/Documents/CUDA/build-EBR_lib-Desktop/');
% addpath('/home/jhagemann/Documents/CUDA/EBR3/build-EBR_lib_santa-Desktop/');
addpath('/home/jhagemann/Documents/CUDA/BeamReconstruction/Matlab/');
addpath(genpath(working_dir))
addpath(genpath(TB_path))
addpath('/home/jhagemann/Documents/CUDA/mmp_sim');

p = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH'); %clears path
p = [p ':/usr/local/cuda/lib'] % new path
setenv('LD_LIBRARY_PATH', p);
clear p

set(0,'DefaultFigureColor','w')
cmap = flipud(colormap('bone(512)'));
set(0,'DefaultFigureColormap', cmap);

cmd = sprintf(['axis image; colorbar;  set(gcf,''PaperPosition'',[0 0 5 5], ''PaperSize'', [5 5]);',...
    ' set(gca,''xtick'',[],''ytick'',[])']);
set(0,'DefaultImageCreateFcn', cmd)
set(groot,'defaultLineLineWidth',1.5)
clear cmap cmd;

%% Simulation parameters
upscale =  8;
p.width = 512;
p.height = p.width;
p.width2 = 1024;
p.height2 = p.width2;
p.rec_width = p.width2;
p.rec_height = p.width2;

p.lambda = 1.24/8; %nm E = 8 keV
p.pixel_size_det = 700; %nm
p.pixel_size = 50   ;
p.F = 1e-3;
p.z1 = 0.32e9;
p.b_0 = 0.99;
p.b_m = 0.75;
p.b_s = 150;
p.beta = 7; % for kaiser bessel window
p.Amp_valid = logical(true(p.rec_height, p.rec_width));

%% Sample Setup
if(1) %dicty sketch
    sample = prepare_probe('dicty_sketch.png', 'dicty_sketch.png',  0, 1, 1, 1, p, 0);
    sample = abs(sample) .* exp(1i.*(angle(conj(sample))));
    p.supp = load('support_dicty_ff.mat');
    p.supp = p.supp.supp.I_supp;
    p.supp = imresize(p.supp, (p.height/256), 'box');
    %     p.supp = pad_to_size(p.supp, p.rec_height, p.rec_width);
end

if(0) % use binary bitmap
    
    sample = binary_bitmap(10, 10, -1, upscale, p.rec_height, p.rec_width);
    p.supp = ones(10*upscale, 10*upscale);
    p.supp(1:floor(10*upscale/2), 1:ceil(10*upscale/2)) = 0;
end

p.supp = pad_to_size(p.supp, p.rec_height, p.rec_width, 'zero');

figure()
imagesc(angle(mid(sample,p)));

figure()
imagesc( mid(p.supp, p));
%% prep data
% holography
prop = PropagatorGPU(p.F, p.F, p.width2, p.height2,1);
id_holo = abs(prop.propTF(sample)).^2;
id_holo = id_holo ./ sum(id_holo(:));
figure; imagesc(id_holo)
% CDI
id_FT = abs(DFT(sample)).^2;
id_FT = id_FT ./ sum(id_FT(:));
figure; imagesc(log10(abs(id_FT)));

%%
% p.num_photons = [1:10:90 100:100:3000];
%  p.num_photons = [1:0.5:20 21:45]; %holo
p.num_photons = [1 25:25:900 1000:500:2500]; %cdi

% p.num_photons = [200];
skip = 1;
% p.num_photons = 300;
repetitons_per_noise = 100;
resolution = zeros(numel(p.num_photons), 2, repetitons_per_noise);
iterations = 200;
cut_off = 15;
gauss = [0];

%%
%  pool = parpool(8);
%% run calculation
tic
for ii = 1:skip:numel(p.num_photons)
    photons = p.num_photons(ii)
    %% holo
    if(0)
        parfor jj = 1:repetitons_per_noise
%         for jj = 1:repetitons_per_noise
            holo = imnoise( (id_holo*photons.* numel(id_holo))*1e-12, 'poisson')*1e12;
            holo = holo ./ sum(holo(:)) .* numel(holo);
            
            holo_rec = RAAR_nf(sqrt(double(holo)), ones(p.rec_height),...
                                iterations, p.F, p);
            
            tmp = angle(holo_rec) .* p.supp;
            
            tmp(tmp < -0.5) = -1;
            tmp(tmp >= -0.5) = 0;
            
            tmp2 = (mid(angle(sample), p)) - (mid(tmp, p));
            
            resolution(ii, 1, jj) = sum(sum(abs(tmp2)));
            
            if(0)                
                figure(1);
                imagesc(mid(tmp, p)); title(sprintf('holo %g photons Δ=%g', photons, resolution(ii, 1, jj)))
                drawnow;
            
%             figure(2); 
%             imagesc(tmp2)
            end
        end
        %         clear x0
    end
    %% CDI
    if(1)
        parfor jj = 1:repetitons_per_noise
%         for jj = 1:repetitons_per_noise
            FT = imnoise((id_FT * photons .* numel(id_FT))*1e-12, 'poisson')*1e12;
            FT(FT == 0) = 2*eps;
            FT = FT ./ sum(FT(:)) .* numel(FT);
            guess = exp(1i .*rand(p.rec_height));
%              guess(~p.supp) = exp(1i*0);
            %             guess = exp(1i .* angle(guess));
%             % psi(~supp) = 1
%             guess(angle(guess) >
%             guess(angle(guess) <= -0.5) = exp(-1i); -0.5) = 1;
            [CDI_rec] = RAAR_ff(sqrt(double(FT)), guess,...
                iterations, p);
            CDI_rec = conj(CDI_rec);
%             [err, aligned_sample] = dftregistration(fft2(angle(CDI_rec)),...
%                 fft2(angle(sample)), 1);
            
%             aligned_sample = (ifft2(aligned_sample));
            
            tmp = angle(CDI_rec) .* p.supp;
            
            tmp(tmp < -0.5) = -1;
            tmp(tmp >= -0.5) = 0;
            
            tmp2 = mid((angle(sample)), p) -...
                                    mid(tmp, p);
            
            resolution(ii, 2, jj) = sum(sum(abs(tmp2)));

            if(0)
                figure(2);
                imagesc(mid(tmp, p));
                title(sprintf('cdi %g photons Δ=%g', photons, resolution(ii, 2, jj)))
                %
                %             figure;
                %             imagesc(tmp2)
            end
            
        end
        %         clear x0 CDI_rec
    end
    
    
end
toc
char(datetime, 'yyyy_MM_dd_''T''HH_mm_ss')

%%
resolution_per_noise = mean(resolution, 3);
std_of_resolution = std(resolution, 0, 3);
%%
figure
hold on
errorbar((p.num_photons), resolution_per_noise(:,1),std_of_resolution(:,1), 'DisplayName', 'holo')
errorbar((p.num_photons), resolution_per_noise(:,2),std_of_resolution(:,2), 'DisplayName', 'CDI')
hold off
% semilogx(p.num_photons, resolution_per_noise(:,1), p.num_photons, resolution_per_noise(:,2))
% hold on
% semilogx(p.num_photons, resolution(:,1), 'DisplayName', 'holo')
% semilogx(p.num_photons, resolution(:,2), 'DisplayName', 'CDI')
% hold off
xlabel('log10(Photons per pixel)')
ylabel('Resolution')
legend show
% ylim([0 0.6])
fname = sprintf('dose_vs_resolution_F_%g_%g_%g_%s.pdf',...
    p.F, p.height, p.rec_height, char(datetime, 'yyyy_MM_dd_''T''HH_mm_ss') )
export_fig(fname)

% save dicty3.mat
