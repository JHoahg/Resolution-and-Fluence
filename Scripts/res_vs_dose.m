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
p.num_iterations = 100;
p.beta = 7; % for kaiser bessel window
p.Amp_valid = logical(true(p.rec_height, p.rec_width));


%% Sample Setup
% for kk = 1:numel(gauss)
if(1) %dicty sketch
    sample = prepare_probe('dicty_sketch.png', 'dicty_sketch.png',  0, 1, 1, 1, p, 0);
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
%    pool = parpool(8);
%%

% p.num_photons = [1:10:90 100:100:2000 3000:1000:20000];
% p.num_photons = [1:10:90 100:100:2000 3000:1000:20000];
%  p.num_photons = [1:0.5:20 21:45]; %holo
% p.num_photons = [1:0.5:20 21:45 70:25:900 1000:500:2500];
p.num_photons = [1:0.5:10 20:5:100 200:100:2000];

skip = 1;
% p.num_photons = 300;
repetitons_per_noise = 50;
resolution = zeros(numel(p.num_photons), 2, repetitons_per_noise, 'distributed');
% reconstructions = zeros( numel(p.num_photons), repetitons_per_noise, p.width * p.height, 'distributed');
iterations = 200;
cut_off = 15;
gauss = [0];
%% run calculation
tic
for ii = 1:skip:numel(p.num_photons)
    photons = p.num_photons(ii)
    %% holo
    if(1)
        parfor jj = 1:repetitons_per_noise
%             for jj = 1:repetitons_per_noise
            %             holo =  makenoisy(id_holo, photons);
            holo = imnoise( (id_holo*photons.* numel(id_holo))*1e-12, 'poisson')*1e12;
            holo(holo == 0) = 2*eps;
            holo = holo ./ sum(holo(:)) .* numel(holo);
            
            holo_rec = RAAR_nf(sqrt(double(holo)), ones(p.rec_height), iterations, p.F, p);
            
%             [err] = dftregistration(fft2(angle(holo_rec)),...
%                 fft2(angle(sample)), 10); err
                        
            
%             figure(1);
%             imagesc(mid(angle(holo_rec), p)); title(sprintf('holo %g photons', photons))
%                         
%             figure(4);
%             imagesc(abs(mid(angle(holo_rec) - angle(sample), p)));
            %
            %             frc_holo_to_ori = FSC(angle((sample)),...
            %                 angle((holo_rec)), p);
            %
            %             [x0, y0, iout,jout ] = intersections(frc_holo_to_ori.nu(cut_off:end),...
            %                 abs(frc_holo_to_ori.frc(cut_off:end)),...
            %                 frc_holo_to_ori.nu(cut_off:end),...
            %                 abs(frc_holo_to_ori.T_hbit(cut_off:end)));
            %
            %             if(isempty(x0) ~= 1)
            %                 resolution(ii, 1, jj) = x0(1);
            %             else
            tmp =sum(sum( abs( angle(mid(holo_rec, p)) - angle(mid(sample, p))).^2));
            resolution(ii, 1, jj) = tmp;
%             tmp
            %             end
            
            %             reconstructions(ii, jj, :) = tmp(:);
            %             resolution(ii, 1, jj)
        end
        %         clear x0
    end
    %% CDI
    if(1)  
        ang_sample = angle(sample);
        parfor jj = 1:repetitons_per_noise
%         for jj = 1:repetitons_per_noise
            %             FT = makenoisy(id_FT, photons);
            FT = imnoise((id_FT * photons .* numel(id_FT))*1e-12, 'poisson')*1e12;
            FT(FT == 0) = 2*eps;
            FT = FT ./ sum(FT(:)) .* numel(FT);
            
            [CDI_rec] = RAAR_ff(sqrt(FT), exp(1i .* 0.2 .*rand(p.rec_height)),...
                iterations, p);
%             CDI_rec = conj(CDI_rec);
            [err] = dftregistration(fft2(angle(CDI_rec)),...
                fft2(ang_sample), 10);
            aligned_sample = circshift_sp(ang_sample, [err(3), err(4)]);
%                                     err
%             aligned_sample = abs(ifft2(aligned_sample));
            %             frc_CDI_to_ori = FSC(angle((CDI_rec)), aligned_sample, p);
            
            %             frc_CDI_to_ori = FSC(angle((CDI_rec)), angle(sample), p);
%             %
%                         figure(2);
%                         imagesc(angle(mid(CDI_rec, p)));
%                         title(sprintf('CDI %g photons', photons))
%             %
%                         figure(3);
%                         imagesc(abs(mid(angle(CDI_rec) - (aligned_sample), p)));
%             %
            
            %             [x0, y0, iout,jout ] = intersections(frc_CDI_to_ori.nu(cut_off:end),...
            %                 abs(frc_CDI_to_ori.frc(cut_off:end)),...
            %                 frc_CDI_to_ori.nu(cut_off:end),...
            %                 abs(frc_CDI_to_ori.T_hbit(cut_off:end)));
            %
            %             if(isempty(x0) ~= 1)
            %                 resolution(ii, 2, jj) = x0(1);
            %             else
            %                 resolution(ii, 2, jj) = 0.5;
            %             end
            %             resolution(ii, 2, jj)
            tmp =sum(sum( abs(angle(mid(CDI_rec, p)) - (mid(aligned_sample, p))).^2));
            resolution(ii, 2, jj) = tmp;
%             tmp
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
%     resolution_per_fwhm{kk} = resolution;        n
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
