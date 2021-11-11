function savedir = print_figures(varargin)
%PRINT_ Summary of this function goes here
%   Detailed explanation goes here

fmt = 'yyyy-MM-dd''T''HH:mm:ss';
load('987b1d18-156e-4597-9811-365a17360dcd.mat', 'fmtext');

p=inputParser;
addRequired(p, 'figs')
addParameter(p, 'folder', tempdir);
addParameter(p, 'base_plot_name', 'plots');
addParameter(p, 'idstr', '');
addParameter(p, 'vstr', datetime('now', 'Format', fmt));
addParameter(p, 'filefmts', {'-dpdf', '-dsvg', '-dpng'});
parse(p, varargin{:});
unpackStruct(p.Results);

fname = sprintf('%s_%s_%s', idstr, base_plot_name, vstr);
fname = strip(fname, '_');

savedir = fullfile(folder, fname);
disp(savedir);
if ~exist(savedir, 'dir')
   mkdir(savedir);
end


pdfnames = cell(numel(figs),1);

for i = 1:numel(figs)
    fig = figs(i);
    if ~isempty(fig.Name)
        filename = fig.Name;
    else
        filename = 'fig';
    end
    filename = sprintf('%s_%s_%s', idstr, base_plot_name, filename);
    
    for j = 1:numel(filefmts)
        temp = nextname(fullfile(savedir, filename), '00', fmtext(filefmts{j}), true);
        print(fig, temp, filefmts{j});    
        fprintf('%s.....saved successfully\n', temp);
        
        if strcmp(filefmts{j}, '-dpdf')
            pdfnames{i} = temp;
        end
    end
    
    
    
end

append_pdfs(fullfile(savedir, [fname '.pdf']), pdfnames{:});
open(fullfile(savedir, [fname '.pdf']));

end

function unpackStruct (structure)
    fn = fieldnames(structure);
    
    for i = 1:numel(fn)
        fni = string(fn(i));
        field = structure.(fni);
        assignin('caller', fni, field);
    end
end



