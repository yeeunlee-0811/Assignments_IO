function [fnvalue,P] = profitFOC(price, theta, inp)
%fnvalue set to be zero to get optimal price

prodcharac = inp.prodcharac;
R = inp.R;
Krc = inp.Krc;
Ntheta1 = inp.Ntheta1_blp;
nu = inp.nu;
marketID = inp.marketID;
PlanID = inp.PlanID;
nM = inp.nM;
nP = inp.nP;
jm = inp.jm;
ablp = inp.ablp;
mcblp = inp.mcblp;

xrc = prodcharac;

theta1 = theta(1:Ntheta1);
sigma = theta(Ntheta1+1:end);

% Recalculate delta
% % To get counterfactual share, we re-calculate delta first
% % Take into account that whatever mc is firms take it as 0.25 less
% % than actual mc
cntf_price = price(PlanID);

X = [ones(jm,1), prodcharac,cntf_price];
delta = X*theta1;

% Recalculate shares
mu = 0; 
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



avgp_uniq = zeros(nP,1);
for p = 1:nP
    ip = PlanID ==p;
    avgp_uniq(p,1) = mean(cntf_price(ip,1));
end

avgshare_uniq = zeros(nP,1);
for p = 1:nP
    ip = PlanID ==p;
    avgshare_uniq(p,1) = mean(P(ip,1));
end

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

elamat_blp = mean(ela_mat,3,'omitnan');
dsdp = diag(elamat_blp);

fnvalue = avgp_uniq + avgshare_uniq./dsdp - (mcblp-(ones(nP,1).*0.25));

