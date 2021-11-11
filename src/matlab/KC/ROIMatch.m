classdef ROIMatch < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nam1 char = ''
        nam2 char = ''
        K1  {mustBeNumeric, mustBeNonnegative, mustBeInteger} = 0
        K2  {mustBeNumeric, mustBeNonnegative, mustBeInteger} = 0
        binC = []
        matched_pairs   (:,2) {mustBeNumeric, mustBePositive, mustBeInteger} = []
        matched_ROIs1   (:,1) {mustBeNumeric, mustBePositive, mustBeInteger} = []
        matched_ROIs2   (:,1) {mustBeNumeric, mustBePositive, mustBeInteger} = []
        nonmatched_ROIs1(:,1) {mustBeNumeric, mustBePositive, mustBeInteger} = []
        nonmatched_ROIs2(:,1) {mustBeNumeric, mustBePositive, mustBeInteger} = []
        unique_match = true
    end
    
    methods
        function obj = ROIMatch(nam1, nam2, varargin)

            p = inputParser;
            p.addRequired('nam1', @(x) isempty(x) || ischar(x));
            p.addRequired('nam2', @(x) isempty(x) || ischar(x));
            p.addParameter('K1', []);
            p.addParameter('K2', []);
            p.parse(nam1, nam2, varargin{:});
            params = p.Results;
            % Copy into properties
            param_fields = fieldnames(params);
            for fn=param_fields(:)'
               obj.(fn{1}) = params.(fn{1}); 
            end            
        end
        
        function init(obj,K1,K2)
            if nargin>1
               obj.K1 = K1;
               obj.K2 = K2;
            end
            
            if obj.K1>0 && obj.K2>0
                obj.binC = false(obj.K1,obj.K2);
            end
            obj.update();
        end
        
        function setBinC(obj, binC)
            obj.binC = binC;
            obj.update();
        end
        
        function match = findMatch(obj, cids1, cids2)
            if isempty(cids1) && ~isempty(cids2)
                binc = obj.binC(:,cids2)';
            elseif isempty(cids2) && ~isempty(cids1)
                binc = obj.binC(cids1,:);
            end
                
            match = rowfun(@find, binc);

        end
        
        function addMatch(obj, varargin)
            if nargin==3
               cids1 =  varargin{1};
               cids2 =  varargin{2};
            elseif nargin==2
                cids1 =  varargin{1}(1);
                cids2 =  varargin{1}(2);
            end
            assert(length(cids1)==length(cids2), 'list of ROIs not the same size')

            ncid = length(cids1);
            
            %clear previous matches if unique_match is true
            if obj.unique_match
               obj.binC(cids1,:) = 0; 
               obj.binC(:,cids2) = 0;
            end
            
            % set new matches
            for i = 1:ncid
               obj.binC(cids1(i), cids2(i)) = 1; 
            end
            obj.update();
        end
        
        function deleteMatch(obj, varargin)
            if nargin==3
               cids1 =  varargin{1};
               cids2 =  varargin{2};
            elseif nargin==2
                cids1 =  varargin{1}(1);
                cids2 =  varargin{1}(2);
            end
            
            if ~isempty(cids1) && ~isempty(cids2)
                assert(length(cids1)==length(cids2), 'list of ROIs not the same size')
                ncids = length(cids1);
                for i = 1:ncids
                   obj.binC(cids1(i), cids2(i)) = 0; 
                end
            elseif ~isempty(cids1)
                obj.binC(cids1,:) = 0;
            elseif ~isempty(cids2)
                obj.binC(:,cids2) = 0;
            end
            obj.update();
        end
        
        function update(obj)
            [K1,K2] = size(obj.binC);
            obj.K1 = K1;
            obj.K2 = K2;
            [r,c] = find(obj.binC);
            obj.matched_pairs = [r c];
            obj.matched_ROIs1 = unique(r);
            obj.matched_ROIs2 = unique(c);
            obj.nonmatched_ROIs1 = setdiff(1:obj.K1, obj.matched_ROIs1);
            obj.nonmatched_ROIs2 = setdiff(1:obj.K2, obj.matched_ROIs2);
        end
        
    end
end

