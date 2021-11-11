function [movie_map, panels, panel_tables] = get_stimulus_info(filepath)
%GET_STIMULUS_INFO Summary of this function goes here
%   Detailed explanation goes here

hallem_set = hallem.load_hallem_set();

%% load txt file, split into odor and movie-stimulus sections

txt = string(fileread(filepath));
txt = splitlines(txt);
txt = strip(txt);
txt(txt=="")=[];
last_panel_line = find(txt=="______________channelA______________");


%% get panel information
panels = hallem.parse_panels(txt, hallem_set);

panel_tables = containers.Map();
keySet = panels.keys();

for i = 1:panels.Count
    hallem_ind = panels(keySet{i}).hallem_inds;
    pin = panels(keySet{i}).pins;
    nam = panels(keySet{i}).nam;
    
    panel = hallem_set.T(hallem_ind,:);
    panel.odor = panel.Properties.RowNames;
    panel.pin = pin;
    panel.Properties.Description = nam;
    panel_tables(nam) = panel;
end
%% get movie information

movie_info = hallem.parse_stimuli(txt);
%channelA = containers.Map();

for i = 1:length(movie_info)
    nam = movie_info(i).nams;
    movie_info(i).T = panel_tables(movie_info(i).panel);
end
%%
movie_map = containers.Map();
for i = 1:numel(movie_info)
    temp = struct();
    temp.nam = movie_info(i).nams;
    temp.panel = movie_info(i).panel;
    temp.channelA = movie_info(i).channelA;
    temp.T = movie_info(i).T;
    
    movie_map(movie_info(i).nams) = temp;
   
end

%end

