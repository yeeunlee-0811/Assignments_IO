% this function gives deltas and choice probabilities P
function [delta, P] = getdelta(sigma,inp)

jm = inp.jm;
nM = inp.nM;
logS = inp.logS;
marketID = inp.marketID;
nu = inp.nu;
prodcharac = inp.prodcharac;
premium = inp.premium;
logit = inp.logit;
Krc = inp.Krc;
R = inp.R;

x = [prodcharac, premium];
xrc = prodcharac;

tol = 1e-5;
dst = 1;
delta = 0.5*ones(jm,1); % arbitrarily chosen starting values for delta

if logit
    while dst > tol
        U = delta;
        eU = exp(U); %jm-by-1
        denom = zeros(nM,1);
        for m = 1:nM
            ix = marketID == m;
            denom(m,1) = 1 + sum(eU(ix,:),1);
        end
        denom = denom(marketID,:);
        P = eU./denom;
        delta_new = delta + logS - log(P);
        dst = max(abs(delta - delta_new));
        delta = delta_new;
    end    
else 
    mu = 0; % util part related to sigmas
    
    for k=1:Krc
        % mu: jm-by-R
        mu = mu + sigma(k)*nu(marketID,:,k).*repmat(xrc(:,k),1,R);
    end
    
    while dst > tol
        U = repmat(delta,1,R) + mu;
        eU = exp(U); %jm-by-R
        denom = zeros(nM,R);
        for m = 1:nM
            ix = marketID == m;
            denom(m,:) = 1 + sum(eU(ix,:),1);
        end
        denom = denom(marketID,:);
        P = mean(eU./denom,2);
        delta_new = delta + logS - log(P);
        dst = max(abs(delta - delta_new));
        delta = delta_new;
    end
end


end




