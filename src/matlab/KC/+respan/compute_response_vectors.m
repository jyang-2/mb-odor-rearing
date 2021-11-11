function response = compute_response_vectors(DFF, responseOptions, ti)
%UNTITLED2 Summary of this function goes here
% Inputs:
%       DFF : [K, T] dF/F traces, one row per neuron/component 
%   responseOptions (optional)
%                 std_thr: threshold for stdev around baseline
%                     thr: threshold for change in dF/F from baseline
%                   max_k: # of maximum values that are averaged for peak value

%        ti : trial/timing information, struct must contain the fields
%               baseline_start
%               baseline_end
%               peakresp_start
%               peakresp_end
%
% Outputs:
%      response
%       .baseline_STD
%       .baseline_MED
%       .peak_response
%       .peakAmp
%       .bin_response
%       .peak_ind
% 
%       .responders
%       .nonresponders
%       .pfo_responders



if ~exist('responseOptions', 'var')
    responseOptions.std_thr = 3;
    responseOptions.thr = 1;
    responseOptions.max_k = 1;
end
% responseOptions.str = sprintf(['std thr = %0.2f\n',...
%                                    'thr = %0.2f\n',...
%                                  'max k = %d'], responseOptions.std_thr, responseOptions.thr, responseOptions.max_k);

responseOptions.str = sprintf('std thr = %0.2f; thr = %0.2f; max k = %d', responseOptions.std_thr, responseOptions.thr, responseOptions.max_k);
responseOptions.str = string(responseOptions.str);
K = size(DFF,1);
    
ti.num_stim = length(ti.peakresp_start);

baseline_STD = NaN(K, ti.num_stim);
baseline_MED = NaN(K, ti.num_stim);
peakMean = NaN(K, ti.num_stim);
peak_response = NaN(K, ti.num_stim);
peak_ind = NaN(K, ti.num_stim);

for it = 1:ti.num_stim
    sdev=std(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),1,2);
    baseline_STD(:,it)=sdev;

    med=median(DFF(:,ti.baseline_start(it):ti.baseline_end(it)),2);
    baseline_MED(:,it) = med;

    [peak_response(:,it), peak_ind(:,it)]=max(DFF(:,ti.peakresp_start(it):ti.peakresp_end(it)),[],2);
    peakMean(:,it) = mean(maxk(DFF(:,ti.peakresp_start(it):ti.peakresp_end(it)), responseOptions.max_k,2),2);
end
peakAmp = peakMean - baseline_MED;
%bin_response = (peakAmp>(baseline_MED+responseOptions.std_thr*baseline_STD)) & (peakAmp>responseOptions.thr);
%bin_response = (peakAmp>(responseOptions.std_thr*baseline_STD)) & (peakAmp>responseOptions.thr);
bin_response = (peakAmp>(responseOptions.std_thr*baseline_STD)) & (peakAmp>responseOptions.thr);
peak_ind = peak_ind + ti.peakresp_start(1:ti.num_stim);

response.options = responseOptions;
response.baseline_STD = baseline_STD;
response.baseline_MED = baseline_MED;
response.peak_response = peak_response;
response.peak_ind=peak_ind;

response.peakAmp = peakAmp;
response.bin_response = bin_response;

%%
responders = (sum(response.bin_response,2)>0);
nonresponders = find(sum(response.bin_response,2)==0);

% pfo_ind = find(cellfun(@(x) strcmp(x, 'pfo'), ti.pin_odors));
% if ~isempty(pfo_ind)
%     pfo_bin_response = response.bin_response(:,ti.si==pfo_ind);
%     pfo_responders = any(pfo_bin_response,2);
% else
%     pfo_responders = false(size(response.bin_response,1),1);
% end
% 
% response.responders = responders;
% response.nonresponders = nonresponders;
% response.pfo_responders = pfo_responders;
%%

%response.bin_corr_resp = corrcoef(bin_response(responders & ~pfo_responders,:), 'Rows', 'pairwise');
%response.bin_corr_all = corrcoef(bin_response(~pfo_responders,:), 'Rows', 'pairwise');
%response.bin_corr_resp = corrcoef(bin_response(responders & ~pfo_responders,:));
%response.bin_corr_all = corrcoef(bin_response(~pfo_responders,:));
%response.bin_corr = response.bin_corr_all;

%response.peakAmp_corr_resp = corrcoef(peakAmp(responders & ~pfo_responders,:), 'Rows', 'pairwise');
%response.peakAmp_corr_all = corrcoef(peakAmp(~pfo_responders,:), 'Rows', 'pairwise');
%response.peakAmp_corr_resp = corrcoef(peakAmp(responders & ~pfo_responders,:));
%response.peakAmp_corr_all = corrcoef(peakAmp(~pfo_responders,:));
%response.peakAmp_corr = response.peakAmp_corr_all;

%[out,idx] = sort(ti.si);            % idx = indices used to sort columns by odor

% response.sorted.bin_corr_resp = corrcoef(bin_response(responders & ~pfo_responders,idx), 'Rows', 'pairwise');
% response.sorted.bin_corr_all = corrcoef(bin_response(~pfo_responders,idx), 'Rows', 'pairwise');
% response.sorted.bin_corr = response.sorted.bin_corr_all;
% 
% response.sorted.peakAmp_corr_resp = corrcoef(peakAmp(responders & ~pfo_responders,idx), 'Rows', 'pairwise');
% response.sorted.peakAmp_corr_all = corrcoef(peakAmp(~pfo_responders,idx), 'Rows', 'pairwise');
% response.sorted.peakAmp_corr = response.peakAmp_corr_all;

% response.sorted.bin_corr_resp = corrcoef(bin_response(responders & ~pfo_responders,idx));
% response.sorted.bin_corr_all = corrcoef(bin_response(~pfo_responders,idx));
% response.sorted.bin_corr = response.sorted.bin_corr_all;
% 
% response.sorted.peakAmp_corr_resp = corrcoef(peakAmp(responders & ~pfo_responders,idx));
% response.sorted.peakAmp_corr_all = corrcoef(peakAmp(~pfo_responders,idx));
% response.sorted.peakAmp_corr = response.peakAmp_corr_all;

%response.sorted.hind_order = ti.panel_table.hallem_ind(out)';
% bin_corr_all
% bin_corr_resp
% peakAmp_corr_all
% peakAmp_corr_resp

end

