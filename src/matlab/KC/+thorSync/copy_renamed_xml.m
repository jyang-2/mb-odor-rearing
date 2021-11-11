function copy_renamed_xml(varargin)
%COPY_RENAMED_THORSYNC Copies and renames all thorSync .xml files from
%  Syncdata folders to an output directory

p = inputParser;
addOptional(p, 'input_folder', []);
addOptional(p, 'output_folder', []);
parse(p, varargin{:});

if isempty(p.Results.input_folder) || isempty(p.Results.output_folder)
   input_folder = uigetdir('../../..'); 
   output_folder = uigetdir('../../..');
end

listing = dir(fullfile(input_folder, '**/ThorRealTimeDataSettings.xml'));

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

for i=1:length(listing)
    filepath = fullfile(listing(i).folder, listing(i).name);
    folders = regexp(filepath, filesep, 'split');
    XmlName = folders(end-1);
    disp([XmlName{1} '.xml'])
    outputBaseFileName = [XmlName{1} '.xml'];
    outputFullFileName = fullfile(output_folder, outputBaseFileName);
    copyfile(filepath, outputFullFileName);
end


end

