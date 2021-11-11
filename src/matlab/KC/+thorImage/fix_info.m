function newThorImageExperiment = fix_info(ThorImageExperiment)

vars = fields(ThorImageExperiment);
bad_fields = {'Date', 'Camera', 'ZStage', 'Timelapse' 'LSM', 'Streaming'};

newThorImageExperiment = ThorImageExperiment;

for i = 1:numel(vars)
    
    fieldname = vars{i};
    if ismember(fieldname, bad_fields) && iscell(ThorImageExperiment.(fieldname)) && numel(ThorImageExperiment.(fieldname))>1
        newThorImageExperiment.(fieldname) = ThorImageExperiment.(fieldname){1};
    end
end
