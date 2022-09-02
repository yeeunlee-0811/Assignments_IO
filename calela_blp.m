function [ela_mat_fin] = calela_blp(sigma,inp)

prodcharac = inp.prodcharac;
R = inp.R;
Krc = inp.Krc;
nu = inp.nu;
marketID = inp.marketID;
PlanID = inp.PlanID;
nM = inp.nM;
nP = inp.nP;
ablp = inp.ablp;

xrc = prodcharac;
delta = getdelta(sigma,inp);

mu = 0; % util part related to sigmas
for k=1:Krc
    % mu: jm-by-R
    mu = mu + sigma(k)*nu(marketID,:,k).*repmat(xrc(:,k),1,R);
end

U = repmat(delta,1,R) + mu;
eU = exp(U); %jm-by-R
denom = zeros(nM,R);
for m = 1:nM
    ix = marketID == m;
    denom(m,:) = 1 + sum(eU(ix,:),1);
end
denom = denom(marketID,:);
P = mean(eU./denom,2);

ela_mat = zeros(nP,nP,nM);
for i = 1:nP
    for j = 1:nP
        for l = 1:nM
            imi = marketID == l & PlanID == i;
            imj = marketID == l & PlanID == j;
            if any(imi)
                if any(imj)
                    if i == j
                        ela_mat(i,j,l) = ablp.*(1-P(imi)).*P(imi);
                    else
                        ela_mat(i,j,l) = -ablp.*P(imj).*P(imi);
                    end
                else
                    ela_mat(i,j,l)=NaN;
                end
            else
                ela_mat(i,j,l)=NaN;
            end
            
        end
    end
end

ela_mat_fin = mean(ela_mat,3,'omitnan');