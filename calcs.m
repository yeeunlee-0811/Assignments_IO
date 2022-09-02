function sum_avgeU = calcs(price, theta,inp)

jm = inp.jm;
nM = inp.nM;
nP = inp.nP;
prodcharac = inp.prodcharac;
premium = inp.premium;
nu = inp.nu;
Krc = inp.Krc;
PlanID = inp.PlanID;
marketID = inp.marketID;
R = inp.R;

price_jm = price(PlanID);

xrc = prodcharac;
X = [ones(jm,1),prodcharac, price_jm];
Ntheta1 = size(X,2);
theta1 = theta(1:Ntheta1);
sigma = theta(Ntheta1:end);

mu = 0; 
for k=1:Krc
    % mu: jm-by-R
    mu = mu + sigma(k)*nu(marketID,:,k).*repmat(xrc(:,k),1,R);
end

delta = X*theta1;

U = repmat(delta,1,R) + mu;
eU_jm = exp(U); %jm-by-R
eU = mean(eU_jm,2); %jm-by-1

avgeU = zeros(nP,1);
for p = 1:nP
    ip = PlanID == p;
    avgeU(p,1) = mean(eU(ip));
end

sum_avgeU = sum(avgeU);



