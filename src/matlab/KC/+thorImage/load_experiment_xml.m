function ThorImageExperiment = load_experiment_xml(xml_filepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
imInfo = xml2struct(xml_filepath);
ThorImageExperiment = imInfo.ThorImageExperiment;
ThorImageExperiment = thorImage.fix_info(ThorImageExperiment);

end

