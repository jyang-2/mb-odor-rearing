function intvl = get_interval_vec(ti)

% check that interval info across trials is consistent
trial_len = unique(ti.trial_end - ti.trial_start + 1);
stim_on = unique(ti.stim_on - ti.trial_start+1);
stim_off = unique(ti.stim_off - ti.trial_start+1);
baseline_start = unique(ti.baseline_start - ti.trial_start+1);
baseline_end = unique(ti.baseline_end - ti.trial_start+1);
peakresp_start = unique(ti.peakresp_start - ti.trial_start+1);
peakresp_end = unique(ti.peakresp_end - ti.trial_start+1);
resp_start = unique(ti.resp_start - ti.trial_start+1);
resp_end = unique(ti.resp_end - ti.trial_start+1);

vars = {trial_len, stim_on, stim_off, baseline_start, baseline_end,...
    peakresp_start, peakresp_end, resp_start, resp_end};
chk = cellfun(@numel, vars);
assert(all(chk), 'Intervals are not consistent across different trials.')

% make logical vectors for indexing 
base_vec = false(1,trial_len);

baseline = base_vec;
baseline(baseline_start:baseline_end) = true;

stim = base_vec;
stim(stim_on:stim_off)=true;

peakresp = base_vec;
peakresp(peakresp_start:peakresp_end) = true;

resp = base_vec;
resp(resp_start:resp_end) = true;

intvl.baseline = baseline;
intvl.stim = stim;
intvl.peakresp = peakresp;
intvl.resp = resp;

    
end

