function single_expt_metadata = load_single_expt_metadata(folder)
%LOAD_SINGLE_EXPT_METADATA Parses metadata for single experiment from 
%                           file 'single_expt_metadata.yaml' in folder
%  
%
% EXAMPLE FILE: 
%-------------------------------------------------------------
%         date: 2021-03-23
%         fly: 1
%         microscope: resonance-galvo
%         thorimagename: movie
%         thorsyncname:
%         olfactometer:
%            pin:  [ 37,  38,  39,  40,  41,  45,  46,  47,  48,  49]
%            hi:   [111,  74,  74,  74,  74, 111,  73,  73,  73,  73]
%            odor:
%            conc: [  0,  -6,  -5,  -4,  -4,   0,  -6,  -5,  -4,  -3]
%         stimulus:
%            conc:     [ 0, 0, 0,-6,-6,-6,-5,-5,-5,-4,-4,-4,-3,-3,-3,-6,-6,-6,-5,-5,-5,-4,-4,-4,-3,-3,-3]
%            channelA: [45,45,45,45,46,46,47,45,47,45,48,48,49,45,49,45,46,46,47,47,45,45,48,48,45,49,49]
%            channelB: [37,37,37,38,38,37,37,39,39,40,37,40,41,41,37,38,37,38,37,39,39,40,37,40,41,37,41]
   
   
fname = fullfile(folder, 'single_expt_metadata.yaml');
fname = char(fname);
single_expt_metadata = ReadYaml(fname);
single_expt_metadata.date = datestr(single_expt_metadata.date);

Hallem = hallem.HallemSet();

stimulus = single_expt_metadata.stimulus;
    stimulus.channelA = cell2mat(stimulus.channelA);
    stimulus.channelB = cell2mat(stimulus.channelB);

    olfactometer = single_expt_metadata.olfactometer;
    olfactometer.pin = cell2mat(olfactometer.pin);
    olfactometer.hi = cell2mat(olfactometer.hi);
    olfactometer.odor = Hallem.odor_abbrev(olfactometer.hi);
    olfactometer.conc = cell2mat(olfactometer.conc);

single_expt_metadata.stimulus = stimulus;
single_expt_metadata.olfactometer = olfactometer;

end

