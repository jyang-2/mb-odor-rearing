 function [data,ThorImageExperiment] = load_experiment_raw(folder)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

xml_filepath = fullfile(folder, 'Experiment.xml');
raw_filepath = fullfile(folder, 'Image_001_001.raw');

imInfo = xml2struct(xml_filepath);
ThorImageExperiment = imInfo.ThorImageExperiment;
ThorImageExperiment = thorImage.fix_info(ThorImageExperiment);


frames = str2double(ThorImageExperiment.Streaming.Attributes.frames);
pixelX = str2double(ThorImageExperiment.LSM.Attributes.pixelX);
pixelY = str2double(ThorImageExperiment.LSM.Attributes.pixelY);
FOV = [pixelY pixelX frames];

    fid = fopen(raw_filepath, 'r');
    if fid == -1
      error('Cannot open file: %s', raw_filepath);
    end
    data = uint16(fread(fid, 'uint16'));
    fclose(fid);
    data = reshape(data, FOV);
    data = pagetranspose(data);

end