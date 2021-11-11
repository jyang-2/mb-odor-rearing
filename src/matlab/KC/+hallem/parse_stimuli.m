function stim_info = parse_stimuli(txt)
last_panel_line = find(txt=="______________channelA______________");
stim_txt = txt(last_panel_line+1:end);
stim_info = struct();

    for i = 1:numel(stim_txt)

       str = stim_txt(i);
       channelA = extractBetween(str, "{", "}");
       channelA = regexp(channelA, '\d*', 'match');
       channelA = str2double(channelA);
       %disp(channelA)
       stim_info(i).channelA = channelA;

       str = eraseBetween(str, "{", "}", 'Boundaries', 'inclusive');
       sstr = strip(split(str,{',', ';', ':', '+'}));
       ind = cellfun(@isempty,sstr, 'UniformOutput', false);

       sstr = sstr(~cell2mat(ind));
%       disp(sstr)

       nams = [];
       for j = 1:numel(sstr)
           exp1 = 'fly\s*(?<fly>\d+)';
           exp2 = 'panel\s*(?<panel_nam>\w+)';
           [~,tok, nonmatch] = regexpi(sstr(j),'(fly|panel)\s*(\w+)', 'match', 'tokens', 'split');

           if ~isempty(tok)      
               stim_info(i).(tok{1}(1)) = tok{1}(2);
           else
               nams = [nams, nonmatch];
           end
       end
       stim_info(i).nams = nams;
    end
    
    for i = 1:length(stim_info)
       stim_info(i).panel = lower( stim_info(i).panel);
    end
end