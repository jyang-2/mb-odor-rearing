function cnm = extract_trial_df_f(CNM, num_trials, fcn)
%EXTRACT_TRIAL_DF_F extract df/f trial by trial
%   Detailed explanation goes here

cnm.T = CNM.T/num_trials;
cnm.C = reshape(full(CNM.C)', [], num_trials, size(CNM.C,1));  % time x trials x cells
cnm.C = permute(cnm.C, [3 1 2]);    % cells x time x trial

cnm.f = reshape(full(CNM.f)', [], num_trials, size(CNM.f,1));  % time x trials x cells
cnm.f = permute(cnm.f, [3 1 2]);    % cells x time x trial
%f1 = reshape(CNM.f(1,:)', [], num_trials); % [time x trial]
%f2 = reshape(CNM.f(2,:)', [], num_trials);
%cnm.f = cat(3, f1, f2); % time x trial x component
%cnm.f = permute(cnm.f, [3 1 2]); %component x time x trial

cnm.R = reshape(full(CNM.R)', [], num_trials, size(CNM.R,1));
cnm.R = permute(cnm.R, [3 1 2]);


cnm.C_df = [];
cnm.F0 = [];
cnm.sC_df = [];
for i = 1:num_trials
    C = cnm.C(:,:,i);
    f = cnm.f(:,:,i);
    %f = cnm.f(1:2,:,i);
    R = cnm.R(:,:,i);
    [C_df,F0] = detrend_df_f(CNM.A,CNM.b,...
        C,f,R,CNM.options);
    sC_df = fcn(C_df);
    %sC_df = smoothdata(C_df, 2, 'loess', 20);
    cnm.C_df = [cnm.C_df C_df];
    cnm.sC_df = [cnm.sC_df sC_df];
    cnm.F0 = [cnm.F0 F0];
end
end

