% Usage:
% 
% savefig(figurehandle,filename)
% savefig(filename)
%
% Saves the current or specified figure with the filename as .pdf, .png and
% .fig file

function save_fig(figurehandle,filename)

	if nargin==1
		filename = figurehandle;
		figurehandle = gcf;
	end
    
	tmp = version;
    if str2double(tmp(1:3))>8.2
        savefig([filename '.fig']);
        export_fig(figurehandle,[filename '.pdf']);
        export_fig(figurehandle,[filename '.png']);
    else
        saveas(figurehandle,filename,'fig');
        export_fig(figurehandle,[filename '.pdf']);
        export_fig(figurehandle,[filename '.png']);   
    end
end

