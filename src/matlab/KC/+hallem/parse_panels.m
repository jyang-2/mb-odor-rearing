function all_panels = parse_panels(txt, hallem)
    last_panel_line = find(txt=="______________channelA______________");
    if isempty(last_panel_line)
       last_panel_line = length(txt)+1; 
    end
    panel_lines = find(contains(txt, "odor panel", 'IgnoreCase', true));
    if ~isempty(last_panel_line)
       panel_lines = panel_lines(panel_lines<last_panel_line); 
    end
    n = numel(panel_lines);
    panel_intervals = [panel_lines; last_panel_line];
    
    panel_names = string(missing);
    panel_conc = [];
    expression = '(\d+)e-?(\d+)';

    % get panel names, concentrations if listed
    for i = 1:n
        str = lower(txt(panel_lines(i)));
        
        [tokens,match] = regexp(str, 'odor panel (\w*)', 'tokens', 'match');
        panel_names(i) = tokens{1};


        [tokens,match] = regexpi(str, expression, 'tokens', 'match');
        if ~isempty(match)
           panel_conc(i)=match;  %#ok<AGROW,
        end
    end

    
    hallem_odors = string(hallem.T.Properties.RowNames);
    hallem_abbrev = string(hallem.T.abbrev);

    all_panels = containers.Map();

    for i = 1:n
        str = txt(panel_intervals(i):(panel_intervals(i+1)-1));
        expression = repmat("\d+\s*:", size(str,1), 1);
        [startIndex,endIndex] = regexp(str,expression);

        str = str(~cellfun(@isempty, startIndex));
        startIndex = cell2mat(startIndex);
        endIndex = cell2mat(endIndex);

        if ~isempty(startIndex)
            pins = nan(size(startIndex));
            odors = strings(size(startIndex));
            abbrevs = strings(size(startIndex));
            hallem_inds = nan(size(startIndex));
            

            for j = 1:numel(pins)
               
               pins(j) = str2double(str{j,1}(startIndex(j):endIndex(j)-1)); 

               substr = strip(str{j,1}((endIndex(j)+1):end));
               %hind = find(cellfun(@(x) contains(substr,x, 'IgnoreCase', true), hallem_odors));
               hind = find(cellfun(@(x) strcmpi(substr, x), hallem_odors));

               if isempty(hind)
                   hind = find(cellfun(@(x) (contains(substr,x, 'IgnoreCase', true) .* (x~="")), hallem_abbrev));
               end           
               if numel(hind)==1
                   hallem_inds(j) = hind;
                   odors(j) = hallem_odors(hind);
                   abbrevs(j) = hallem_abbrev(hind);
               end
                
              %disp([num2str(pins(j)) ':    ' substr ', ' num2str(hind)]);

            end
            
            
            
            panel.pins = pins;
            panel.hallem_inds = hallem_inds;
            panel.odors = odors;
            panel.abbrevs = abbrevs;
            panel.nam = panel_names(i);
            panel.conc = panel_conc;

        end
        all_panels(panel_names(i)) = panel;
    end

end