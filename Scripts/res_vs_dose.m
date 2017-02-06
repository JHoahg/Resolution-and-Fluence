%builds upon dose_cdi_vs_holo.m
%JH 20160725

%% Simulation parameters
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
p.num_iterations = 200;
p.beta = 7; % for kaiser bessel window
p.Amp_valid = logical(true(p.rec_height, p.rec_width));
p.oversample = 1;
% control output & calculations
p.direct_propagation = 1;
p.do_phase_reconstruction = 0;
p.show_results = 0;
p.measure_resolution_FRC = 1;
p.measure_delta = 0;
%% Sample Setup
% for kk = 1:numel(gauss)
if(0) %dicty sketch
    sample = prepare_probe('dicty_sketch.png', 'dicty_sketch.png',  0, 0.1, 1, 1, p, 0);
    sample = abs(sample) .* exp(1i.*(angle(conj(sample))));
    p.supp = load('support_dicty_ff.mat');
    p.supp = p.supp.supp.I_supp;
    p.supp = imresize(p.supp, (p.height/256), 'box');
    %     p.supp = pad_to_size(p.supp, p.rec_height, p.rec_width);
end

if(0) % use binary bitmap
    upscale =  10;
    sample = binary_bitmap(10, 10, -1, upscale, p.rec_height, p.rec_width);
    p.supp = ones(10*upscale, 10*upscale);
    p.supp(1:floor(10*upscale/2), 1:ceil(10*upscale/2)) = 0;
%     p.width = 110;
%     p.height = p.width;
end


if(1) %mandrill dürer
    sample = prepare_probe('mandrill.png', 'durer.png',  -0.4, 0.4, 0.8, 1.2, p, 0);
    sample = abs(sample) .* exp(1i.*(angle(conj(sample))));
    p.supp = ones(size(sample));
    
end
p.supp = pad_to_size(p.supp, p.rec_height, p.rec_width, 'zero');

figure()
imagesc(angle(mid(sample,p)));

figure()
imagesc( mid(p.supp, p));
%% prep data
% holography
prop = PropagatorGPU(p.F, p.F, p.width2, p.height2,1);
id_holo = prop.propTF(sample);
id_holo_phases = angle(id_holo);
id_holo = abs(id_holo).^2;
id_holo = id_holo ./ sum(id_holo(:));
figure; imagesc(id_holo)

% CDI
id_FT = DFT(gpuArray(sample));
id_FT_phases = (angle(id_FT));
id_FT = abs(id_FT).^2;
id_FT = id_FT ./ sum(id_FT(:));
figure; imagesc(log10(abs(id_FT)));
%%
   pool = parpool(8);
%% set up scan parameters

% p.num_photons = [1:10:90 100:100:2000 3000:1000:20000];
% p.num_photons = [1:0.5:20 21:45 70:25:900 1000:500:2500];
% p.num_photons = [1:0.5:10 20:5:100 200:100:2000];
p.num_photons = [10.^[0.1:0.1:4.5]];

skip = 1;

p.repetitons_per_noise = 2;
% resolution = zeros(numel(p.num_photons), 2, p.repetitons_per_noise, 'distributed');
resolution = zeros(numel(p.num_photons), 2, p.repetitons_per_noise);
% reconstructions = zeros( numel(p.num_photons), repetitons_per_noise, p.width * p.height, 'distributed');
iterations = 200;
cut_off = 15;
gauss = [0];
%% run calculation
tic
for ii = 1:skip:numel(p.num_photons)
    photons = p.num_photons(ii)
    ang_sample = angle(mid(sample,p));
    %% holo
    if(1)
%         if(p.direct_propagation)
%             inv_prop = PropagatorGPU(-p.F, -p.F, p.width2, p.height2,1);
%         end
%                 
%         parfor jj = 1:p.repetitons_per_noise
         for jj = 1:p.repetitons_per_noise
            holo = imnoise( (id_holo*photons.* numel(id_holo))*1e-12, 'poisson')*1e12;
            holo(holo == 0) = 2*eps;
            keyboard;
            holo = holo ./ sum(holo(:)) .* numel(holo);

            if(p.direct_propagation)
                inv_prop = PropagatorGPU(-p.F, -p.F, p.width2, p.height2,1);
                holo = sqrt(holo) .* exp(1i .* id_holo_phases);
                holo_rec = gather(inv_prop.propTF(holo));
            end
            
            if(p.do_phase_reconstruction)
                holo_rec = RAAR_nf(sqrt(double(holo)), ones(p.rec_height), iterations, p.F, p);
            end
            
            if(p.show_results)
                figure(1);
                imagesc(mid(angle(holo_rec), p)); title(sprintf('holo %g photons', photons))
                
                figure(4);
                imagesc(abs(mid(angle(holo_rec) - angle(sample), p)));
                drawnow;
            end
            
            if(p.measure_resolution_FRC)
                
                frc_holo_to_ori = FSC(ang_sample,...
                    angle(mid(holo_rec,p)), p);
                
                [x0, y0, iout, jout ] = intersections(frc_holo_to_ori.nu(cut_off:end),...
                    abs(frc_holo_to_ori.frc(cut_off:end)),...
                    frc_holo_to_ori.nu(cut_off:end),...
                    abs(frc_holo_to_ori.T_hbit(cut_off:end)));
                
                if(isempty(x0) ~= 1)
                    resolution(ii, 1, jj) = x0(1);
                else
                    resolution(ii, 1, jj) = 0.5;
                end
            end
            
            if(p.measure_delta)
                tmp = sum(sum( abs( angle(mid(holo_rec, p)) - mid(ang_sample, p) ).^2));
                resolution(ii, 1, jj) = tmp;
            end
        end
    end
    
    %% CDI
    if(1)
        parfor jj = 1:p.repetitons_per_noise
%         for jj = 1:p.repetitons_per_noise
            
            FT = imnoise((id_FT * photons .* numel(id_FT))*1e-12, 'poisson')*1e12;
            FT(FT == 0) = 2*eps;
            FT = FT ./ sum(FT(:)) .* numel(FT);
            
            if(p.direct_propagation)
                FT = sqrt(FT) .* exp(1i .* id_FT_phases);
                CDI_rec = gather(IDFT(FT));
            end
            
            if(p.do_phase_reconstruction)
                [CDI_rec] = RAAR_ff(sqrt(FT), exp(1i .* 0.2 .*rand(p.rec_height)),...
                    iterations, p);
            end
            
            if(p.measure_resolution_FRC)
                
                [err] = dftregistration(fft2(angle(mid(CDI_rec,p))),...
                    fft2(ang_sample), 10);
                aligned_sample = circshift_sp(ang_sample, [err(3), err(4)]);
                
                frc_CDI_to_ori = FSC(angle(mid(CDI_rec,p)), aligned_sample, p);
                                
                [x0, y0, iout,jout ] = intersections(frc_CDI_to_ori.nu(cut_off:end),...
                    abs(frc_CDI_to_ori.frc(cut_off:end)),...
                    frc_CDI_to_ori.nu(cut_off:end),...
                    abs(frc_CDI_to_ori.T_hbit(cut_off:end)));
                
                if(isempty(x0) ~= 1)
                    resolution(ii, 2, jj) = x0(1);
                else
                    resolution(ii, 2, jj) = 0.5;
                end
            end
            
            if(p.measure_delta)
                
                tmp =sum(sum( abs(angle(mid(CDI_rec, p)) - (mid(aligned_sample, p))).^2));
                resolution(ii, 2, jj) = tmp;
            end
            
            if(p.show_results)
                figure(2);
                imagesc(angle(mid(CDI_rec, p)));
                title(sprintf('CDI %g photons', photons))
                drawnow;
                if(p.measure_resolution_FRC)
                    figure(3);
                    imagesc(abs(mid(angle(CDI_rec),p) - mid(aligned_sample, p)));
                end
            end
            
        end
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
errorbar(log10(p.num_photons), resolution_per_noise(:,1),std_of_resolution(:,1), 'DisplayName', 'holo')
errorbar(log10(p.num_photons), resolution_per_noise(:,2),std_of_resolution(:,2), 'DisplayName', 'CDI')
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

% save bitmap_new_poisson.mat
