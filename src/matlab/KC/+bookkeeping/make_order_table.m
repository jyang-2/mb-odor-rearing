function order_table = make_order_table(a, b)
%MAKE_ORDER_TABLE 
% Returns a table for keeping track of odor pairs by hallem index
%
% INPUTS
%    make_order_table(a)
%             a : hpair_ind (index for odor pair
%
%    make_order_table(hind1, hind2)
%         hind1 : hallem index for 1st odor in odor pair
%         hind2 : hallem index for 2nd odor in odor pair
%
% OUTPUTS

hallem_set = hallem.load_hallem_set();

switch nargin
    case 1
        hpair_ind = a;
        
    case 2
        fcn = @(x,y) (hallem_set.pairs.A(x,y));
        hpair_ind = arrayfun(fcn, a, b);
end

hpair_ind = unique(hpair_ind);
pair_cats = categorical(hpair_ind, ...
                        1:height(hallem_set.odors.pairs),...
                        hallem_set.odors.pairs.pair_str, ....
                        'Ordinal', true);

order_table = table();
order_table.hpair_ind = hpair_ind;
order_table.hind1 = hallem_set.odors.pairs.hind1(order_table.hpair_ind);
order_table.hind2 = hallem_set.odors.pairs.hind2(order_table.hpair_ind);
order_table.pair_cats = pair_cats;

% pred = rowfun(fcn, order_table, 'InputVariables', {'hind1', 'hind2'},...
%                                  'OutputVariableNames', {'model_corr', 'model_cos', 'model_jaccard'});
% order_table.model_corr = pred.model_corr;
% order_table.model_cos = pred.model_cos;
% order_table.model_jaccard = pred.model_jaccard;

end

