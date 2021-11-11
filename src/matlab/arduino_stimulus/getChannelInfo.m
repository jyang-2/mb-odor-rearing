function stimulus = getChannelInfo(channelA, olfactometer)
%GETCHANNELINFO Takes stimulus pin list and olfactometer struct, and returns 
%                 struct with stimulus info, 

% inputs:
%   - channelA : [1 x nStim] Row vector of arduino pins used to deliver stimuli
%   - olfactometer : [struct]
%         - pin : [1 x nPin] row vector with arduino pins, used for
%                  indexing the other variables and channelA lookup
%         -  hi : [1 x nPin] row vector with hallem index
%         - odor: [1 x nPin] string array with name of odors hooked up
%         - conc: [1 x nPin] row vector with odor concentrations
%         -    T: [nStim x 4] table with variables {'pin', 'hi', 'odor', 'conc'}
% outputs:
%   - stimulus: [struct]
%       - pin : [1 x nStim] arduino pins (=channelA from input)
%       - hi : [1 x nStim] row vector w/ hallem index of delivered odors
%       - odor : [1 x nStim] string array, row vector of delivered odors
%       - conc : [1 x nStim] row vector with odor concentrations
%

arguments
    channelA (1,:)
    olfactometer struct
end
[~,idx] = ismember(channelA, olfactometer.pin);

pin = channelA;
hi = olfactometer.hi(idx);
odor = olfactometer.odor(idx);
conc = olfactometer.conc(idx);

stimulus.pin = pin;
stimulus.hi = hi;
stimulus.odor = odor;
stimulus.conc = conc;
stimulus.T = table(pin', hi', odor', conc', 'VariableNames', {'pin', 'hi', 'odor', 'conc'});
end

