%% generate hallem pair lookup table

count = 0;
A = zeros(111);

for i = 1:111
    for j = i:111    
        count = count+1;
        A(i,j) = count;
    end
end

A = (A+A') - eye(size(A,1)).*diag(A);

[hind2, hind1] = meshgrid(1:111, 1:111);

hpair_ind = A(:);
hind1 = hind1(:);
hind2 = hind2(:);

hallem_pair_lookup = table(hind1, hind2, hpair_ind);
writetable(hallem_pair_lookup, 'hallem_pair_lookup.csv');
%%
ordered_pairs = struct();
count = 0;
for i = 1:111
    for j = i:111    
        count = count + 1;
        ordered_pairs(count).hind1 = i;
        ordered_pairs(count).hind2 = j;
        ordered_pairs(count).o1 = hallem_set.odors.odor_list(i);
        ordered_pairs(count).o2 = hallem_set.odors.odor_list(j);
        ordered_pairs(count).abbrev1 = hallem_set.odors.abbrev(i);
        ordered_pairs(count).abbrev2 = hallem_set.odors.abbrev(j);
    end
end
%%
hallem_set = hallem.load_hallem_set();
%%
load('hallem_pair_mat.mat')
%%
o1 = hallem_set.odors.pairs.o1;
o2 = hallem_set.odors.pairs.o2;
pair_str = strcat(o1, ", ", o2);
hallem_set.odors.pairs.pair_str = pair_str;

