function modelKC = get_model_kcs()
    hs = readtable('hs100.csv');
    modelKC.jaccard = zeros(110);
    modelKC.smc = zeros(110);
    modelKC.corr = zeros(110);
    modelKC.cos = zeros(110);
    
    for i = 1:height(hs)
        modelKC.jaccard(hs.o1(i), hs.o2(i)) = hs.jaccard(i);
        modelKC.smc(hs.o1(i), hs.o2(i)) = hs.smc(i);
        modelKC.corr(hs.o1(i), hs.o2(i)) = hs.corr(i);
        modelKC.cos(hs.o1(i), hs.o2(i)) = hs.cos(i);
    end

end

