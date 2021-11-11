classdef HallemSet < handle
    %HALLEMSET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T
        chem
        mat_filename = 'hallem-8ea3408c-de8c-42a6-a609-8bf4c3770264'
        cmap_filename = 'favorite_colormaps.mat';
        odors
        pairs
        receptor
        cmaps
        valid_only = true;
    end
    
    methods
        function obj = HallemSet(mat_filename)
            if nargin>0
                obj.mat_filename = mat_filename;
            end
            temp = load(obj.mat_filename);
            cmaps = load(obj.cmap_filename);
            obj.T = temp.T;
            obj.chem = temp.chem;
            obj.odors = temp.odors;
            obj.pairs = temp.pairs;
            obj.receptor = temp.receptor;
            obj.cmaps = cmaps;
        end
        
        function good = good_odor(obj, hi)
            %GOOD_ODORS returns logical array for valid odors
            if isnumeric(hi) && all(round(hi)==hi, 'all')
                good = ismember(hi, 1:111);
            elseif ischar(hi) || isstring(hi) || iscellstr(hi)
                hi = string(hi);
                good_full = ismember(hi, obj.odors.odor_list);
                good_abbrev = ismember(hi, obj.odors.abbrev(obj.odors.abbrev~=""));
                
                if any(good_full)
                    good = good_full;
                elseif any(good_abbrev)
                    good = good_abbrev;
                else
                    good = false(size(hi));
                end
            end
        end
        
        function names = odor_name(obj,hi)
            %ODOR_NAME returns the full odor name of index or abbreviation
            %   Detailed explanation goes here
            if obj.valid_only
                assert(all(obj.good_odor(hi),'all'), 'Invalid odors found');
            end

            odor_names = obj.odors.odor_list;
            if isnumeric(hi) && all(round(hi)==hi, 'all')
                names = strings(size(hi));
                [good, idx] = ismember(hi, 1:111);
                names(good) = odor_names(idx(good));
           elseif ischar(hi) || isstring(hi) || iscellstr(hi)
               hi = string(hi);
               names = strings(size(hi));
               [good, idx] = ismember(hi, obj.odors.abbrev);
               names(good) = odor_names(idx(good));
            end
        end
        
        function index = odor_index(obj, hi)
            if obj.valid_only
                assert(all(obj.good_odor(hi),'all'), 'Invalid odors found');
            end

            if ~isstring(hi)
                hi = string(hi);
            end
            index = zeros(size(hi));
            [good_abbrev, idx_abbrev] = ismember(hi, obj.odors.abbrev);
            [good_name, idx_name] = ismember(hi, obj.odors.odor_list);
            
            index(good_abbrev) = idx_abbrev(good_abbrev);
            index(good_name & ~good_abbrev) = idx_name(good_name & ~good_abbrev);
            
        end
        
        function abbrev = odor_abbrev(obj, hi)
            if obj.valid_only
                assert(all(obj.good_odor(hi),'all'), 'Invalid odors found');
            end
            %ODOR_ABBREV returns the abbreviation of an odor's index or name
            odor_abbrev = obj.odors.abbrev;
            if isnumeric(hi) && all(round(hi)==hi, 'all')
                abbrev = strings(size(hi));
                [good, idx] = ismember(hi, 1:111);
                abbrev(good) = odor_abbrev(idx(good));
           elseif ischar(hi) || isstring(hi) || iscellstr(hi)
               hi = string(hi);
               abbrev = strings(size(hi));
               [good, idx] = ismember(hi, odor_abbrev);
               abbrev(good) = odor_names(idx(good));
            end
        end
        
        function hpi = odor_pair_index(obj, hi1, hi2)
            if obj.valid_only
                assert(all(obj.good_odor(hi1),'all'), 'Invalid odors found');
                assert(all(obj.good_odor(hi2),'all'), 'Invalid odors found');
            end
            assert(all(size(hi1) == size(hi2),'all'), 'Hallem indices are different sizes');
            
            A = obj.pairs.A;
            hpi = zeros(size(hi1));
            for i = 1:numel(hi1)
                if all(ismember([hi1(i) hi2(i)], 1:111))
                    hpi(i) =  A(hi1(i), hi2(i));
                end
            end
        end
        
        function [hi1, hi2] = hpi_odor_index(obj, hpi)
            Tpair = obj.pairs.T;
            
           [good_pairs, idx] = ismember(hpi, Tpair.hpair_ind);
            if obj.valid_only
               assert(all(good_pairs, 'all'), 'Invalid odor pair values found');
            end
           hi1 = zeros(size(hpi));
           hi2 = zeros(size(hpi));
           
           hinds = [Tpair.hind1(idx(good_pairs)) Tpair.hind2(idx(good_pairs))];
           
           hi1(good_pairs) = min(hinds,[],2);
           hi2(good_pairs) = max(hinds,[],2);
        end
    end
    
    methods(Static)
       function iPN_homo = PNtransform(ORNc)
           %UNTITLED Summary of this function goes here
           %   Detailed explanation goes here
           % PN = rmax * (ORN^1.5) / (ORN^1.5 + s^1.5 + sigma^1.5)
           % s = m_inp*EAG/190
           
           % parameters for transformation
           rmax = 165;
           n = 1.5;
           sigma = 12;
           m_inp = 10.63;
           
           % zero all negative firing rates
           if istable(ORNc)
              ORN =  ORNc{:,:};
           elseif ismatrix(ORNc)
               ORN = ORNc;
           end
           ORN(ORN<0)=0;
           
           % compute the EAG based on raw Hallem ORN values
           EAG = sum(ORN,2);
           EAG(EAG<0) = 0;
           
           % compute suppression factor
           s = m_inp*EAG/190;
           
           % input gain, homogeneous, set at VM7 levels from integration study
           if istable(ORNc)
               iPN_homo = ORNc;
               iPN_homo{:,:} = rmax*(ORN.^n)./((ORN.^n)+sigma^n+s.^n);
           elseif ismatrix(ORNc)
               iPN_homo = rmax*(ORN.^n)./((ORN.^n)+sigma^n+s.^n);
           end
           
           
       end

 
    end
end

