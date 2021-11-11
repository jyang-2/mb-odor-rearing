function olfactometer = cell2olfactometer(arduino_pin_info)
%CELL2OLFACTOMETER Takes a cell array [m x 4] and returns a struct called
% "olfactometer" containing pin, odor, and concentration information.
%
% 
% inputs:
%   - arduino_pin_info : pin data for olfactometer in cell array, with columns
%                        corresponding to {pin, hallem index, odor, concentration}
%                         
%                         it's okay if the hallem index column or the odor
%                         column is empty
% outputs:
%   - olfactometer : struct with the fields
%         - pin : [1 x nPin] row vector with arduino pins
%         -  hi : [1 x nPin] row vector with hallem index
%         - odor: [1 x nPin] string array with name of odors hooked up
%         - conc: [1 x nPin] row vector with odor concentrations
%         -    T: [nStim x 4] table with variables {'pin', 'hi', 'odor', 'conc'}
% 


    arguments
        arduino_pin_info (:,4)
    end
    
    pin = cell2mat(arduino_pin_info(:,1));
    hi = cell2mat(arduino_pin_info(:,2));
    odor = string(arduino_pin_info(:,3));
    conc = cell2mat(arduino_pin_info(:,4));
    T = table(pin, hi, odor, conc, 'VariableNames', {'pin', 'hi', 'odor', 'conc'});
    
    olfactometer = struct();
    olfactometer.pin = pin';
    olfactometer.hi = hi';
    olfactometer.odor = odor';
    olfactometer.conc = conc';
    olfactometer.T = T;
    
    
 %  input odor pin data in cell array, with each row being 
 %  {pin, hallem index, concentration}
end

