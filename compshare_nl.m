function fitshare = compshare_nl(beta,price,inp)

jm = inp.jm;
marketID = inp.marketID;
PlanID = inp.PlanID;
nest = inp.nest;
prodcharac = inp.prodcharac;
lninNS = inp.lninNS;

Xnl = [prodcharac, lninNS, price];
inp.Xnl=Xnl;

sigma = beta(5,1);

X = [ones(jm,1), Xnl];
%deltajm: jm-by-1 vector of delta
xi = residualFun_nl(beta,inp);
deltajm = X*beta +xi;

% exp_delta: jm-by-1 vector that goes in to Dgm expression
exp_delta = exp(deltajm/(1-sigma));

Dgm = zeros(jm,1);
for i = 1:jm
    m = marketID(i);
    n = nest(i);
    imn = marketID == m & nest ==n; %imn gives index of products which are in the same nest&market
    Dgm(i,1) = sum(exp_delta(imn,1));
end

% sum_dgm is a component of denominator which is the same within each
% market
sum_dgm = zeros(jm,1);
for i = 1:jm
    m = marketID(i);
    im1 = marketID ==m & nest == 1;
    im0 = marketID ==m & nest == 0;
    dgm_n1 = mean(Dgm(im1,1)).^(1-sigma);
    if isnan(dgm_n1)
        dgm_n1 = 0;
    end
    dgm_n0 = mean(Dgm(im0,1)).^(1-sigma);
    if isnan(dgm_n0)
        dgm_n0 = 0;
    end
    sum_dgm(i,1) = dgm_n1 + dgm_n0 +1;
end
checksumdgm = [marketID, PlanID,nest, Dgm,sum_dgm];

denom = (Dgm.^sigma).*sum_dgm;

fitshare = exp_delta./denom;


    
    



