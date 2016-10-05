% hio_fresnel performs iterative reconstruction only dependend on the
% fresnel number
%
%	Usage:
%	[psi, result]	= hio_fresnel(im,fresnelx,fresnely,settings)
%	defaultsettings = hio_fresnel
%  
%	im			= intensity 
%	fresnelx	= fresnelnumber for horizontal image frequencies
%	fresnely	= fresnelnumber for vertical image frequencies
%	settings	= structure with other settings
%	  
%	
%	Default settings:
%		beta = 0.2;				% HIO Parameter
%		gamma = 0.2;			% HIO Parameter
%		C = 0;					% maximaler positiver Phasenschub (>pi zum
%		D = 0.01;				% Fehlerkriterium
%		D_break = 1;			% Fehlerkriterium als Abbruchbedingung nutzen
%								deaktivieren)
%	
%		N_it = 500;				% maximale Schrittanzahl
%		N_start = 500000;		% Index ab der beta/delta Kopplung kommt, wenn
%								gr????er als N_it ist es Klaus' Variante
%		ER_steps = 0;			% nach N_start alle ER_steps Schritte au??erhalb
%								des Supports hart auf 1 setzen
%		beta_delta = 0.0019;	% Verh??ltnis von beta zu delta, erst ab N_start
%								benutzt
%		plotlive = 0;			% W??hrend der Iteration alle plotlive Schritte
%								plotten
%		statusint = 25;			% Interval in dem eine Statusmeldung (inkl.
%								vergangener Zeit) kommen soll
%		startsol = []			% startsol = 'holo' f??r holographische
%								Startl??sung oder ein komplexes Feld ??bergeben
%		a0  = 1;				% F??r Zufallsstartl??sungen Parameter mit
%		phi0 = 1;				% denen multipliziert wird
%	
%		save_int = 0;			% Intervall in dem gespeichert wird
%		save_filename = './rec';  % Dateiname zum speichern
%

function [psi, result]=hio_fresnel(im,fresnelx,fresnely,settings)
funstart = tic;
if (nargin < 4)
	settings = struct;
end

%% Default settings

% enforce padding factor
defaults.padfactor = 0;

% rebinning of images to reduce calc. time
defaults.binfactor = 0;

% use GPU if images are small enough
defaults.use_gpu = 0;

% dont use support or unwrapp info
defaults.supp = 0;
defaults.phi_unwrap = 0;
defaults.do_unwrap = false;


% Defaults:
defaults.a0 = 1;
defaults.phi0= 0;

defaults.beta_delta = 0;	
defaults.abs_seq_one = false;
defaults.beta = 0.8;
defaults.gamma = 0.8;
defaults.C = 4;
defaults.D = 0;
defaults.D_break = true;

defaults.plotlive = 100;
defaults.plotonly = 0;
defaults.plotcolorbar = 1;
defaults.N_it = 10000;
defaults.showtime = true;

defaults.statusint = 25;
defaults.save_int = 0;
defaults.save_filename = './rec.mat';

defaults.use_ER_ampl = 0;
defaults.use_ER_phase = 0;
defaults.noamplconstraint = 0;
defaults.upd_only_outside_padarea = false;
defaults.upd_only_outside_padarea_fourierconstraint = false;
defaults.use_realspacesamplingpropagator = 0;
defaults.use_stepbystep = 0;
defaults.purephase_assumption = 1;

defaults.startsol = 'random_klaus';
defaults.method = 'hio';


fields=fieldnames(defaults);

% return default settings
if (nargin == 0)
	psi=defaults;
	return
end

% optionale settings aus struct auslesen
for i=1:size(fields,1)
	thisfield = fields{i};
	if (isfield(settings,thisfield))
		curset.(thisfield) = settings.(thisfield);
	else
		curset.(thisfield) = defaults.(thisfield);
	end
end 

[Ny,Nx] = size(im);
px = abs(1/(Nx*fresnelx));
py = abs(1/(Ny*fresnely));


%% dependend settings
if curset.padfactor==0
	curset.padfactor = max([px,py,1]);
end

if curset.supp == 0
	curset.supp = ones(Ny,Nx);
end

if curset.phi_unwrap==0
	curset.phi_unwrap=zeros(Ny,Nx);
end

% binning
if curset.binfactor>1
    fprintf('Images are binned by a factor of %.2f to improve SNR\n',curset.binfactor);
% keyboard;
    im					= imresize(im,1/curset.binfactor);
    curset.supp			= imresize(curset.supp,1/curset.binfactor);
    curset.phi_unwrap	= imresize(curset.phi_unwrap,1/curset.binfactor);
    
	if ~strcmpi(curset.startsol,'random') && ~strcmpi(curset.startsol,'holo') && ~strcmpi(curset.startsol,'random_klaus')
		curset.startsol = padarray(curset.startsol,[pady padx],'replicate');
		curset.startsol = imresize(curset.startsol,1/curset.binfactor);
	end
    fresnelx = fresnelx * curset.binfactor^2;
    fresnely = fresnely * curset.binfactor^2;
    [Ny, Nx] = size(im);
end


% pad if needed or if the user wanted it
padx = 0;
pady = 0;
if curset.padfactor>1
	fprintf('Images are enlarged by a factor of %.2f to ensure propper propagation \n',curset.padfactor);
	pady = round((ceil(curset.padfactor)*Ny -Ny)/2);
	padx = round((ceil(curset.padfactor)*Nx -Nx)/2);

	im   = padarray(im,[pady padx],'replicate');
	curset.supp = padarray(curset.supp,[pady padx],'replicate');
	curset.phi_unwrap = padarray(curset.phi_unwrap,[pady padx],'replicate');
	
	if ~strcmpi(curset.startsol,'random') && ~strcmpi(curset.startsol,'holo') && ~strcmpi(curset.startsol,'random_klaus')
		curset.startsol = padarray(curset.startsol,[pady padx],'replicate');
	end
	
	% new image size
	[Ny, Nx] = size(im);
	
	% check if GPU is still ok	
	if curset.use_gpu
		if (max(Ny,Nx) < 5000)
			curset.use_gpu = 1;
		else
			fprintf('\t Tried to use GPU but images are to large. \n');
			curset.use_gpu = 0;
		end
	end
end


[Ny, Nx] = size(im);

% Propagator Kerne
if curset.use_realspacesamplingpropagator==1
	
	% coordinate system in real space
	[X, Y] = meshgrid((1:Nx)-floor(Nx/2)-1,(1:Ny)-floor(Ny/2)-1);

	% convolution kernel in real space
	prop_back	= fft2c(exp(-1i*pi*(fresnelx*X.^2+fresnely*Y.^2)))*sqrt(fresnelx*fresnely);
	prop_forward = fft2c(exp(+1i*pi*(fresnelx*X.^2+fresnely*Y.^2)))*sqrt(fresnelx*fresnely);
else
	% Koordinatensystem im Fourierraum 
	dqx = 2*pi/(Nx);
	dqy = 2*pi/(Ny);
	[Qx,Qy] = meshgrid(((1:Nx)-floor(Nx/2)-1)*dqx,((1:Ny)-floor(Ny/2)-1)*dqy);
	
	% Propagatorkern
	kappa_x = -1/2*(Qx.^2);
	kappa_y = -1/2*(Qy.^2);
	prop_back	= exp(-1i*kappa_x/(2*pi*fresnelx)-1i*kappa_y/(2*pi*fresnely)); 
	prop_forward = exp(+1i*kappa_x/(2*pi*fresnelx)+1i*kappa_y/(2*pi*fresnely));
end

% Amplituden
Ampl = sqrt(im);

% initial guess
if (isfield(curset,'startsol'))
	if strcmpi(curset.startsol,'holo')
		% holographische 
		Psi = Ampl;
		psi = ifft2c(fft2c(Psi).*prop_back);
	elseif strcmpi(curset.startsol,'random_klaus')
		% random guess
		psi= (1 + curset.a0*(rand(size(im))-0.5)) .* exp(1i*(curset.phi0*(rand(size(im))-0.5)).*(1-curset.supp)  + 1i*  (curset.phi0 + curset.phi0*(rand(size(im))-0.5)).*curset.supp);
	elseif strcmpi(curset.startsol,'random')
		% random guess with right amplitudes		
		Psi = Ampl.*exp(1i*curset.phi0*(rand(size(im))-0.5));
		psi = ifft2c(fft2c(Psi).*prop_back);
	else
		% benutzerdefinierte 
		psi = curset.startsol;
		disp('Individuelle Startloesung wird benutzt');
	end
end

% Werte aus den Argumenten in settings schreiben
settings.startsol = psi;
err = zeros(curset.N_it,1);

% Variablen in GPU schieben
if curset.use_gpu
	psi = gpuArray(psi);
	curset.supp = gpuArray(curset.supp);
	prop_back = gpuArray(prop_back);
	prop_forward = gpuArray(prop_forward);
	Ampl = gpuArray(Ampl);
	curset.phi_unwrap = gpuArray(curset.phi_unwrap);
	err = gpuArray(err);
end

%% Ab hier geht es eig. erst los

ind = 0; 
psi_upd = angle(psi);

tic 

for ind =1:curset.N_it 
	% In Detektorebene propagiert:
	Psi = ifft2c(fft2c(psi).*prop_forward);
		
	% Diskrepanzprinzip, nur innerhalb des urspr�nglichen Bereichs berechnen
	dev2 = (abs(Psi(pady+1:end-pady,padx+1:end-padx)).^2-Ampl(pady+1:end-pady,padx+1:end-padx).^2).^2;
	d = sqrt(sum(dev2(:))/numel(dev2));
	err(ind) = d;
	
	%% Statusbericht
	if (mod(ind,curset.statusint)==0)
		blub = toc;
		fprintf('Iteration #% 4.f  - Fehler: %01.4f - Zeitdifferenz %03.1f s\n',ind,d,blub);
		tic
	end	
	
	%% Diskrepanzprinzip
	if d > curset.D
		Psi2 = (1-curset.D/d).*Ampl.^2 + curset.D/d*abs(Psi).^2;
		% "Fourier Constraint"

		if (curset.upd_only_outside_padarea_fourierconstraint)
			Psi(pady+1:end-pady,padx+1:end-padx) = sqrt(Psi2(pady+1:end-pady,padx+1:end-padx)) ...
				.*Psi(pady+1:end-pady,padx+1:end-padx)./abs(Psi(pady+1:end-pady,padx+1:end-padx));
		else
			Psi = sqrt(Psi2).*Psi./abs(Psi);
		end
	else
		fprintf ('Diskrepanz-Abbruchbedingung nach %i Iterationen erreicht\n',ind);
		Psi2 = (1-curset.D/d).*Ampl.^2 + curset.D/d*abs(Psi).^2;
		
		% "Fourier Constraint"
		if (curset.upd_only_outside_padarea)
			Psi(pady+1:end-pady,padx+1:end-padx) = sqrt(Psi2(pady+1:end-pady,padx+1:end-padx)) ...
				.*Psi(pady+1:end-pady,padx+1:end-padx)./abs(Psi(pady+1:end-pady,padx+1:end-padx));
		else
			Psi = sqrt(Psi2).*Psi./abs(Psi);
		end
		if (curset.D_break || d<curset.D_break)
			% Speichern und dann abbrechen
			res.psi = gather(psi);
			res.err = gather(err(1:ind));
			res.settings = curset;
			result.phase_unwr = phase_unwr;
			if (curset.save_int > 0)
				filename = sprintf('%s_final.mat',curset.save_filename);
				save(filename,'res');
			end
			
			break; 
		end
	end	  

	%% in Objektebene propagiert
    
    psi_upd_old = psi;
	psi_upd = ifft2c(fft2c(Psi).*prop_back); 

	if strcmpi(curset.method,'hio')
        % phase constraint - modified HIO as in K. Giewekemeyer et. al.
        % (Waveguide with dicty)
        upd_phase_o = (angle(psi) - curset.gamma*angle(psi_upd)).*(1-curset.supp);
        upd_phase_i = min(curset.C,angle(psi_upd)).*curset.supp;
        if curset.use_ER_phase
            upd_phase_o = 0.*(1-curset.supp);
        end
        upd_phase = upd_phase_i + upd_phase_o;

        % unwrap phase if set 
        phase_unwr = upd_phase;
        if curset.do_unwrap
            differ = (curset.phi_unwrap-upd_phase)/(2*pi);
            phase_unwr = upd_phase+2*pi*round(differ);
        end

        if curset.noamplconstraint
            abs_psi = abs(psi_upd);
        elseif (curset.beta_delta>0 || curset.beta==0)
            % amplitude constraint - coupling between phase and amplitude as in 
            % M. Krenkel et. al. (HoloTIE) and outside supp HIO constraint    
            % for beta_delta == 0 and beta == 0 its the classical GS constraint
            abs_psi_o = (abs(psi) - curset.beta*(abs(psi_upd) - 1)).*(1-curset.supp);
            abs_psi_i = exp(curset.beta_delta*phase_unwr).*curset.supp;
        else
            % modified HIO scheme by Klaus Giewkemeyer, for absorption its
            % the classical HIO by Fienup
            abs_psi_o = (abs(psi) - curset.beta*(abs(psi_upd) - 1)).*(1-curset.supp);
            if (curset.purephase_assumption)
                abs_psi_i = (abs(psi) - curset.beta*(abs(psi_upd) - 1)).*curset.supp;
            else
                abs_psi_i = abs(psi_upd).*curset.supp;
            end
        end
        
        if curset.use_ER_ampl
            abs_psi_o = 1.*(1-curset.supp);
        end
        abs_psi = abs_psi_i + abs_psi_o;
        
        % Absorption darf nicht > 1 sein:
        if (curset.abs_seq_one)
            abs_psi(abs_psi(:)>1) = 1;
        end
        
        % put together phase and amplitude
        if (curset.upd_only_outside_padarea)
            psi = ones(size(abs_psi),'like',abs_psi);
            psi(pady+1:end-pady,padx+1:end-padx) = abs_psi(pady+1:end-pady,padx+1:end-padx) ...
                .* exp(1i*phase_unwr(pady+1:end-pady,padx+1:end-padx));
        else
            psi = abs_psi .* exp(1i*phase_unwr);
        end        
    elseif strcmpi(curset.method,'gs')
        % GS type real space constraint
        psi = psi_upd./abs(psi_upd);
        
        % unwrap phase if set 
        phase_unwr = angle(psi);
        if curset.do_unwrap
            differ = (curset.phi_unwrap-upd_phase)/(2*pi);
            phase_unwr = upd_phase+2*pi*round(differ);
        end
    end
    

	
	%% Speichern
	if ((mod(ind,curset.save_int) == 0 || ind == curset.N_it) && curset.save_int > 0)
		% Zwischenergebnisse alle paar Schritte und am Ende abspeichern
		res.psi = gather(psi(pady+1:end-pady,padx+1:end-padx));
		res.err = gather(err(1:ind));
		res.settings = settings;
		
		filename = sprintf('%s_%06.f',curset.save_filename,ind);
		save(filename,'res');
	end
	
	%% Plotten
	if (curset.plotlive > 0 && mod(ind,curset.plotlive)==0)
		plot_it_now();
		if (curset.use_stepbystep)
			pause();
		end
	end
	
	if (ind == curset.N_it)	  
		fprintf ('Maximale Anzahl Iterationen (%i) erreicht.\n',ind); 
	end
end % for

% Ergebnis
psi = gather(psi(pady+1:end-pady,padx+1:end-padx));
result.psi = psi(pady+1:end-pady,padx+1:end-padx);
result.psi_udp = psi_upd_old;
result.err = gather(err(1:ind));
result.phase_unwr = phase_unwr;
result.settings = settings;

plot_it_now();

% Zeit anzeigen
ttoc = toc(funstart);
if (curset.showtime)
	fprintf ('Iterative Rekonstruktion dauerte %g sekunden \n',ttoc);
end


%% helper functions
	function val = fft2c(data)
		val = fftshift(fft2(ifftshift(data)));
	end

	function val = ifft2c(data)
		val = fftshift(ifft2(ifftshift(data)));
	end


	function plot_it_now()
		if curset.plotonly == 0
			% plot all
			subplot(2,3,1)
			imagesc(angle((psi))); title(['angle(psi) - Step:',num2str(ind)]);
			axis equal tight; plotthecolorbar();
			set(gcf,'color','white');
			set(gca,'xtick',[]); set(gca,'ytick',[]);
			
			subplot(2,3,4)
			imagesc(abs((psi))); title(['abs(psi) - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
			set(gcf,'color','white');
			set(gca,'xtick',[]); set(gca,'ytick',[]);
			
			subplot(2,3,2)
			imagesc(angle((psi_upd))); title(['angle(psi)  after constraint  - Step:',num2str(ind)]);
			axis equal tight; plotthecolorbar();
			set(gcf,'color','white');
			set(gca,'xtick',[]); set(gca,'ytick',[]);
			
			subplot(2,3,5)
			imagesc(abs((psi_upd))); title(['abs(psi) after constraint - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
			set(gcf,'color','white');
			set(gca,'xtick',[]); set(gca,'ytick',[]);
			
			subplot(2,3,6)
			imagesc(abs((Psi).^2)); title(['rec. Intensity - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
			set(gcf,'color','white');
			set(gca,'xtick',[]); set(gca,'ytick',[]);
			
			subplot(2,3,3)
			semilogy((err(1:ind))); title(['d - Step:',num2str(ind)]);
			colormap bone(255);
			drawnow;
		else
			switch(curset.plotonly)
				case 1
					imagesc(angle((psi(pady+1:end-pady,padx+1:end-padx)))); title(['angle(psi) after constraint - Step:',num2str(ind)]);
					axis equal tight; plotthecolorbar();
					set(gcf,'color','white');
					set(gca,'xtick',[]); set(gca,'ytick',[]);
				case 2
					imagesc(angle((psi_upd(pady+1:end-pady,padx+1:end-padx)))); title(['angle(psi) - Step:',num2str(ind)]);
					axis equal tight; plotthecolorbar();
					set(gcf,'color','white');
					set(gca,'xtick',[]); set(gca,'ytick',[]);
				case 3
					plot((err(1:ind))); title(['d - Step:',num2str(ind)]);
				case 4
					imagesc(abs((psi(pady+1:end-pady,padx+1:end-padx)))); title(['abs(psi) after constraint - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
					set(gcf,'color','white');
					set(gca,'xtick',[]); set(gca,'ytick',[]);
				case 5
					imagesc(abs((psi_upd(pady+1:end-pady,padx+1:end-padx)))); title(['abs(psi) - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
					set(gcf,'color','white');
					set(gca,'xtick',[]); set(gca,'ytick',[]);
				case 6
					imagesc(abs((Psi(pady+1:end-pady,padx+1:end-padx)).^2)); title(['rec. Intensity - Step:',num2str(ind)]); axis equal tight; plotthecolorbar();
					set(gcf,'color','white');
					set(gca,'xtick',[]); set(gca,'ytick',[]);
			end
			colormap bone(255);
			drawnow;
		end
	end

	function plotthecolorbar()
		if curset.plotcolorbar == 1
			colorbar;
		else
			colorbar off;
		end
	end
end