function tuning = get_tuning(bin_response, peakAmp, hi)
%GET_TUNING tuning = get_tuning(response, hi)
%   returns a structure containing:
%         hi      [1 x ntrials], hallem indices for all trials
%         nhi     # of different odors in hi
%         K       # of cells (from response input struct)
%         uhi     unique hallem indices
%         occ     # of occurrences for each odor in nhi
%
%         ntrials     [K x nhi], # of trials in which cells respond, for each
%                                odor
%         fractrials  [K x nhi], fraction of trials in which cells respond, 
%                                for each odor
%         maxPeakResp [K x nhi], maximum peak amplitude of each responding 
%                                response (if a cell never responds to a 
%                                certain odor, maxPeakResp=0)
%         meanPeakAll [K x nhi], mean peak amplitude of every cell's
%                                response to a given odor 
%
%   function dependencies:
%           uniquecount

%bin_response = response.bin_response;
%peakAmp = response.peakAmp;
binPeakAmp = bin_response.*peakAmp;

[uhi, occ] = uniquecount(hi);
K = size(bin_response,1);
nhi = numel(uhi);

tuning.hi = hi;
tuning.nhi = nhi;
tuning.K = K;
tuning.uhi = uhi;
tuning.occ = occ;

tuning.ntrials = zeros(K,nhi);
tuning.fractrials = zeros(K,nhi);
tuning.maxPeakResp = zeros(K,nhi);
tuning.meanPeakAll = zeros(K,nhi);


for idx = 1:numel(uhi)
    tuning.ntrials(:,idx) = sum(bin_response(:,hi==uhi(idx)),2);   
    tuning.maxPeakResp(:,idx) = max(binPeakAmp(:,hi==uhi(idx)),[],2);
    tuning.meanPeakAll(:,idx) = mean(peakAmp(:,hi==uhi(idx)),2);
end
tuning.fractrials = tuning.ntrials./tuning.occ;

tuning.breadth = sum(tuning.fractrials > 0, 2);
tuning.nodors = sum(tuning.fractrials>0,2);

[~, I] = max(tuning.maxPeakResp, [], 2);
tuning.bestOdor = I;
tuning.bestOdor(tuning.nodors==0)=0;
end

