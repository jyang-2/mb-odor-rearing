classdef TrialTrace
%TRIALTRACE Class for easy extraction of trial traces (select by any combo
%           of cell id, timepoint, and trial).
%
%   trialObj = TrialTrace(C_df, ti)
%           or 
%   trialObj = TrialTrace(C_df_trial) if you already ran trace.trace_to_trial
%
%   trialObj.getTrials('Cells', cids, 'Trials', [1 5 10]) 

    properties
        C_df_trial
        K
        T
        num_trials        
    end
    
    methods
        function obj = TrialTrace(C_df, ti)
            if nargin == 1
               obj.C_df_trial = C_df;
            elseif nargin == 2
                obj.C_df_trial = trace.trace_to_trial(C_df, ti);
            end
            [obj.K,obj.T,obj.num_trials] = size(obj.C_df_trial);
        end
        
        function traces = getTrials(obj,varargin)
            
            p = inputParser;
            p.addParameter('Cells', []);
            p.addParameter('Trials', []);
            p.addParameter('Timepoints', []);
            p.parse(varargin{:});
            params = p.Results;
            
            cids = params.Cells;
            trials = params.Trials;
            ts = params.Timepoints;
            
            % select cells
            traces = obj.C_df_trial;
            if ~isempty(cids)
                traces = traces(cids,:,:);
            end
            
            % select timepoints
            if ~isempty(ts)
                traces = traces(:,ts,:);
            end
            
            % select trials
            if ~isempty(trials)
               traces = traces(:,:,trials); 
            end
        end
    end
end

