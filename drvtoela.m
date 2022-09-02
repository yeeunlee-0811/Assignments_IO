function elamat_fin =  drvtoela(elamat,inp)

jm = inp.jm;
nP = inp.nP;
nM = inp.nM;
PlanID = inp.PlanID;
marketID = inp.marketID;
share = inp.avgshare; % 16x1
price = inp.avgp;% 16x1

elamat_fin = zeros(nP,nP);

for r = 1:nP
    for c = 1:nP
        elamat_fin(r,c) = elamat(r,c)*price(c)/share(r);
    end
end




