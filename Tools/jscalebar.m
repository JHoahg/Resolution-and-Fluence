function len = jscalebar(length, pixelsize, posx, posy, color,sbtext)
% posx and posy in percent of the image!
% example:
% jscalebar(5000, pixelsize, 0.95, 0.95 , 'w');
%SCALEBAR creates a scalebar in the actual figure
%
% Usage: scalebar(length, pixelsize, posx, posy)
    if nargin<3
        posx = 0.90; 
        posy = 0.90;
    end
    if nargin<5
        color = 'k';
	end
	if nargin<6
		sbtext = [num2str(length/1000), ' μm'];
	end
    len = round(length/pixelsize); % 100µm scale bar
    bla = axis;
    posx = bla(2)*posx;
    posy = bla(4)*posy;
    
    line([posx posx-len],[posy posy],'LineWidth',4,'Color',color);
    
    text(posx-0.5*len , posy - bla(2)*0.04, sbtext,...
        'Color',color , 'HorizontalAlignment', 'center')
end

