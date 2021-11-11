%pfo_hex = 
%pfo_pentanol = 
pentanol_conc = [0, -3, -4, -5, -6];
pentanol_pins = containers.Map(pentanol_conc,[nan, nan, nan, nan, nan]);
pentanol_pins(0) = 45;  % pfo
pentanol_pins(-3) = 49;
pentanol_pins(-4) = 48;
pentanol_pins(-5) = 47;
pentanol_pins(-6) = 46;

hexanol_conc = [0, -3, -4, -5, -6];
hexanol_pins = containers.Map(hexanol_conc,[nan, nan, nan, nan, nan]);
hexanol_pins(0) = 37;    % pfo
hexanol_pins(-3) = 41;
hexanol_pins(-4) = 40;
hexanol_pins(-5) = 39;
hexanol_pins(-6) = 38;
%%
conc_list = [0, -3, -4, -5, -6];
pentanol = containers.Map(conc_list, cell(size(conc_list)));
hexanol = containers.Map(conc_list, cell(size(conc_list)));
pfo = [pentanol_pins(0); hexanol_pins(0)];
mix = containers.Map(conc_list, cell(size(conc_list)));

for conc = conc_list
    pentanol(conc)  = [pentanol_pins(conc); hexanol_pins(0)];
    hexanol(conc) = [pentanol_pins(0); hexanol_pins(conc)];
    mix(conc) = [pentanol_pins(conc); hexanol_pins(conc)];
end

%%
conc_block_3 = [pfo, pentanol(-3), hexanol(-3), mix(-3)];
conc_block_4 = [pfo, pentanol(-4), hexanol(-4), mix(-4)];
conc_block_5 = [pfo, pentanol(-5), hexanol(-5), mix(-5)];
conc_block_6 = [pfo, pentanol(-6), hexanol(-6), mix(-6)];
%% SINGLE CONCENTRATION
conc = -3;
conc_block = [pfo, pentanol(conc), hexanol(conc), mix(conc)];
num_blocks = 3;

channelA_blocks = cell(1,num_blocks);
channelB_blocks = cell(1,num_blocks);
for i = 1:num_blocks
   idx = randperm(length(conc_block));
   temp = conc_block(:,idx);
   channelA_blocks{i} = temp(1,:);
   channelB_blocks{i} = temp(2,:);
end

a = cell2mat(channelA_blocks);
b = cell2mat(channelB_blocks);

strA = regexprep( mat2str(a), {'\[', '\]', '\s+'}, {'', '', ','});
strA = ['channelA[] = {' strA '};'];

strB = regexprep( mat2str(b), {'\[', '\]', '\s+'}, {'', '', ','});
strB = ['channelB[] = {' strB '};'];
disp(['-------------------shuffled conc. blocks: ' num2str(conc) '-------------------']) ;
disp(strA);
disp(strB);
%% enter list of concentrations for blocks manually

%conc_list = [-6 -5 -4 -6 -5 -4, -6, -5, -4];
conc_list = [-6,-5,-4,-3];
num_blocks = length(conc_list);

channelA_blocks = cell(1,num_blocks);
channelB_blocks = cell(1,num_blocks);
for i = 1:num_blocks
    conc = conc_list(i);
    
    conc_block = [pfo, pentanol(conc), hexanol(conc), mix(conc)];
    idx = randperm(length(conc_block));
    temp = conc_block(:,idx);
    
    channelA_blocks{i} = temp(1,:);
    channelB_blocks{i} = temp(2,:);
end

a = cell2mat(channelA_blocks);
b = cell2mat(channelB_blocks);

strA = regexprep( mat2str(a), {'\[', '\]', '\s+'}, {'', '', ','});;
strA = ['channelA[] = {' strA '};'];

strB = regexprep( mat2str(b), {'\[', '\]', '\s+'}, {'', '', ','});
strB = ['channelB[] = {' strB '};'];
disp(['-------------------shuffled conc. blocks-------------------']) ;
disp(strA);
disp(strB);
%% pid pins
temp = [pentanol(-6) pentanol(-5) pentanol(-4) pentanol(-3),...
            hexanol(-6) hexanol(-5) hexanol(-4) hexanol(-3),...
            mix(-6) mix(-5) mix(-4) mix(-3)];
a = temp(1,:);
b = temp(2,:);

strA = regexprep( mat2str(a), {'\[', '\]', '\s+'}, {'', '', ','});;
strA = ['channelA[] = {' strA '};'];

strB = regexprep( mat2str(b), {'\[', '\]', '\s+'}, {'', '', ','});
strB = ['channelB[] = {' strB '};'];
disp(['-------------------pid blocks-------------------']) ;
disp(strA);
disp(strB);