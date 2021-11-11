function [new_ti, idx] = ti_frame_corrector(ti)
    new_ti = ti;
    
    idx = zeros(1,length(ti.frameouts.rT));
    for i = 1:length(ti.block_ict)
        idx = idx + (ti.frameouts.rT > ti.block_ict(i) & ti.frameouts.rT < ti.block_fct(i));
    end
    idx = (idx == 0); % "invalid" frames

    if sum(idx) > 0
        disp(['Removing ', num2str(sum(idx)), ' frames (', ...
            num2str(length(ti.frameouts.rT)),' frames to ', ...
            num2str(length(ti.frameouts.rT) - sum(idx)),' frames)']);
        disp(' ');
        disp('Invalid frame time(s):');

        bad_rT = ti.frameouts.rT(idx);
        bad_fT = ti.frameouts.fT(idx);
        bad_T = ti.frameouts.T(idx);
        for i = 1:length(bad_rT)
            disp([num2str(bad_rT(i)), ' to ', num2str(bad_fT(i)), ...
                ' (', num2str(bad_fT(i) - bad_rT(i)),')'])
        end
    end
    new_ti.frameouts.rT = ti.frameouts.rT(~idx);
    new_ti.frameouts.fT = ti.frameouts.fT(~idx);
    new_ti.frameouts.T = ti.frameouts.T(~idx);
    new_ti.frame_times = ti.frame_times(~idx);

    % ------------------------------------------------------------------------------
    % 4. Get odor frame times
    % ------------------------------------------------------------------------------
    % plot(ai.frameouts.Time, ai.Frame_Out,'k'); hold on;
    % plot(ai.frameouts.Time, ai.olfDispPin,'r');

%     idx = zeros(1, length(new_ti.frameouts.rT));
%     for i = 1:length(new_ti.stim_ict) % Use center time
%         idx = idx + (new_ti.frameouts.T > new_ti.stim_ict(i) & new_ti.frameouts.T < new_ti.stim_fct(i));
%     end
%     % idx = (idx == 0); % odor frames
% 
%     new_ti.olf_pulse_rT = new_ti.frameouts.rT(logical(idx));
%     new_ti.olf_pulse_fT = new_ti.frameouts.fT(logical(idx));
%     new_ti.olf_pulse_T = new_ti.frameouts.T(logical(idx));
    
    % Update some info
    new_ti.timepoints = length(new_ti.frameouts.rT);
%    new_ti.ts = length(new_ti.frameouts.rT);
    
    disp(' ');
    disp(['old ti, # of timepoints: ', num2str(length(ti.frame_times))]);
    disp(['new ti, # of timepoints: ', num2str(length(new_ti.frame_times))]);
   
