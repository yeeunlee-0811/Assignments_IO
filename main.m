clear

%% Import data
data0 = importdata('data.xlsx');
data = data0.data;

%% 
marketID = data(:,1);
PlanID = data(:,2);
coverage = data(:,3);
network = data(:,4);
satisfaction = data(:,5);
premium = data(:,6); % price
share = data(:,7);
logS = log(share);

jm = size(data,1); % the mumber of market x product combinations

prodcharac = [coverage, network,satisfaction]; 
npc = size(prodcharac,2); % the number of product characteristics

%% 
uniqMarket = sort(unique(marketID),'ascend'); % total 600 markets
nM = size(uniqMarket,1);
uniqPlan = sort(unique(PlanID),'ascend'); % total 16 plans
nP = size(uniqPlan,1);
%%
% generating dummies for market and products
% % market dummy
Mdummy = zeros(jm,nM); 
for i=1:nM
    m = uniqMarket(i);
    Mdummy(:,i) = marketID == m;
end
Mdummy(:,1)= []; % drop a market for outside good

% % plan dummy
Pdummy = zeros(jm,nP); 
for i=1:nP
    p = uniqPlan(i);
    Pdummy(:,i) = PlanID == p;
end
Pdummy(:,1)= []; % drop a market for outside good

% % the number of plan for each market
nPM = zeros(nM,1);
for i =1:nM
    ip = marketID == i;
    nPM(i,1) = sum(ip);
end

nPMv = zeros(jm,1);
for j = 1:nM
    im = marketID == j;
    nPM0 = nPM(j,1);
    nPMv(im,1) = nPM0;
end

avgp = zeros(jm,1);
for i=1:jm
    p = PlanID(i);
    ip = PlanID == p;
    avgp(i,1) = mean(premium(ip));
end

    
%% Generating instruments - BLP
% BLP instruments: Measure of product's isolation in characteristic space

instr0 = prodcharac; % don't use price to generate instruments 
ninst = size(instr0,2);

Z0 = zeros(jm,ninst,3); 
for i= 1:jm
    market = marketID(i);
    plan = PlanID(i) ; 
    samemarket = marketID == market;
    sameplan = PlanID == plan;
    for j =1:ninst
        tinst = instr0(:,j);
    Z0(i,j,1) = tinst(i);
    Z0(i,j,2) = sum(tinst(samemarket))-tinst(i);
    Z0(i,j,3) = sum(tinst(sameplan))-tinst(i);
    end
end

Z1 = reshape(Z0,jm,ninst*3);






%% Making input structure 
inp.jm = jm;
inp.nM = nM;
inp.nP = nP;
inp.logS = logS;
inp.prodcharac = prodcharac;
inp.marketID = marketID;
inp.PlanID = PlanID;
inp.premium = premium;
inp.share = share;
%% BLP _related input
R = 150; % number of random draws
inp.R = R;
Kbeta = 1+size(prodcharac,2)+1; % number of beta-parameters: constant term + prodcharac + price
Krc = Kbeta-2; % the first Krc prod.chars. have random coeff. Drop constant term & price
inp.Krc = Krc;
% different draws for each year, and each product characteristic.
rng(10) 
nu0 = randn(nM,R,Krc); %Market level draws

nu = zeros(jm,R,Krc);
for i = 1:nM
    im = marketID == i; 
    nui = nu0(i,:,:);
    nu(im,:,:) = ones(sum(im),1).*nui;
end
inp.nu =nu;

%% Logit

inp.logit =1;

theta0 = zeros(1,1);

truedelta= getdelta(theta0,inp);

Z_logit= Pdummy;
inp.Z = Z_logit;

    
IV1 = fitlm(Z_logit,premium);

phat = [ones(jm,1),Z_logit]*IV1.Coefficients.Estimate; 

IV2 = fitlm([prodcharac,phat],truedelta);

%% Export_logit
varname_logit = ["constant", "coverage", "network","satisfaction","premium"]';
logit_estimate = table;
logit_estimate.name = varname_logit;
logit_estimate.coeff = IV2.Coefficients.Estimate;
logit_estimate.se = IV2.Coefficients.SE;
logit_estimate.tStat = IV2.Coefficients.tStat;
logit_estimate.pVlaue = IV2.Coefficients.pValue;

save('logit_estimate.mat','logit_estimate')
save('logit_estimate.txt','logit_estimate')
%% % Logit_Average Elasticiity

% price coefficient of logit 
al = IV2.Coefficients.Estimate(end,1); 
inp.al = al;

elamat_logit = calela_logit(inp);

%% Logit mean marginal cost and markup 

avgp_uniq = zeros(nP,1);
for p = 1:nP
    ip = PlanID ==p;
    avgp_uniq(p,1) = mean(premium(ip,1));
end

avgshare_uniq = zeros(nP,1);
for p = 1:nP
    ip = PlanID ==p;
    avgshare_uniq(p,1) = mean(share(ip,1));
end

inp.avgshare = avgshare_uniq;
inp.avgp = avgp_uniq;
%%
% elasticity is dS/dP*p/S (price is along the row, share is along the
% column)
elamat_fin_l = drvtoela(elamat_logit,inp);
ownpela_l = diag(elamat_logit);
mclogit = avgp_uniq + (avgshare_uniq./ownpela_l);
markup_logit = 100*(avgp_uniq-mclogit)./avgp_uniq; 


%% Nested logit

% Setting groups: 
% Group0 = coverage <=0.8 
% Group1 = coverage > 0.8

nest = [];
in1 = coverage > 0.8;
nest(in1,1) = 1;
inp.nest = nest;

% calculate outside good(not insured) share for each market
% share0: outside good share of each market 
share0 = zeros(jm,1);
for m=1:nM
    im = marketID == m;
    share0(im,:) = 1-sum(share(im,:));
end

% Nshare: Nest share of each market
Nshare = zeros(jm,1);
for i =1:jm
    in = nest(i,1);
    im = marketID(i,1);
    inm = nest == in & marketID == im; %inm gives index of goods which are in the same nest & in the same market
    Nshare(i,1) = sum(share(inm));
end

% lns_s0
lnshare = log(share);
lnshare0 = log(share0);
lns_s0=lnshare-lnshare0;

inp.lnshare= lnshare;
inp.lnshare0=lnshare0;
inp.lns_s0 = lns_s0;

% calculate each plan's share in the nest
inNshare = share./Nshare;
lninNS =log(inNshare);
inp.lninNS = lninNS;

% check 
checkshare = [marketID, nest, share, Nshare, inNshare,share0];




%% nested logit estimation

Z_nlogit = [Pdummy, Nshare];

IV1_nlogit_ns = fitlm(Z_nlogit,lninNS);
IV1_nlogit_p = fitlm(Z_nlogit, premium);

lninNS_hat = [ones(jm,1), Z_nlogit]*IV1_nlogit_ns.Coefficients.Estimate;
phat = [ones(jm,1),Z_nlogit]*IV1_nlogit_p.Coefficients.Estimate;

x_nl = [prodcharac, lninNS_hat,phat];

IV2_nlogit = fitlm(x_nl,lns_s0);
%% Export
varname_nl = ["constant", "coverage", "network","satisfaction","log(share_in_nest)","premium"]';
NLogit_estimate = table;
NLogit_estimate.name = varname_nl;
NLogit_estimate.coeff = IV2_nlogit.Coefficients.Estimate;
NLogit_estimate.se = IV2_nlogit.Coefficients.SE;
NLogit_estimate.tStat = IV2_nlogit.Coefficients.tStat;
NLogit_estimate.pVlaue = IV2_nlogit.Coefficients.pValue;

save('NLogit_estimate.mat','NLogit_estimate')

%% check fitted shares
beta_nl = IV2_nlogit.Coefficients.Estimate;
fitshare = compshare_nl(beta_nl,premium,inp);
dist_max = max(abs(share-fitshare));
dist_min = min(abs(share-fitshare));
dist_mean = mean(abs(share-fitshare));

%% Calculate NL elasticity 
% To calculate NL elasticity, I used function derivative 

elamat_nl = calela_nl(beta_nl,inp);
%%
elamat_fin_nl = drvtoela(elamat_nl,inp);
ownpela_nl = diag(elamat_nl);
mcnl = avgp_uniq + (avgshare_uniq./ownpela_nl);
markup_nl = 100*(avgp_uniq-mcnl)./avgp_uniq; 

%% avgprodcharac
avgprodcharac = zeros(nP,3);
for p = 1:nP
    ip = PlanID == p;
    avgprodcharac(p,:) = mean(prodcharac(p,:),1);
end

%% BLP

% Test "getdelta" function
inp.logit = 0;
sigma =zeros(Krc,1);
testdelta = getdelta(sigma,inp);
%% BLP instruments & Refine

Z_blp0 = [ones(jm,1), Z1, nPMv, avgp] ;

z0_D = Z_blp0(:,1:2);
for i=3:size(Z_blp0,2)
    zi = Z_blp0(:,i);
    reg = fitlm(z0_D,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi);
        end
        z0_D = [z0_D zi];
        % disp(i)
    end
end

Z_blp = z0_D;
inp.Z_blp=Z_blp;
%%

W = (Z_blp'*Z_blp)/jm; % initial weighting matrix
inp.W = W;
inp.conc = true; % "concentrate out" linear parameters
inp.logit = 0;


residual_fn = @(theta) residualFun(theta,inp); 
obj = @(beta) GMM_objective(beta,Z_blp,W,residual_fn);

theta0 = rand(Krc,1); % starting values; randomly drawn from (0,1)

%% Run fminunc
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
theta2_0 = fminunc(obj,theta0,options); % theta2 = sigma, #Krc
[~,theta1_0]=residualFun(theta2_0,inp); % theta1 = parameters in the mean util, # 1+#prodcharac+premium
theta_0 = [theta1_0;theta2_0]; 
inp.Ntheta1_blp = size(theta1_0);
%% Run fminsearch
% options = optimset('Display','iter','MaxFunEvals',20000);
% theta2_0 = fminsearch(obj,theta0,options);

%% Calculate BLP std error
inp.conc =0;
[blpvar] = se_blp(theta_0,inp);
blpse = sqrt(diag(blpvar));
tstat_blp = theta_0./blpse;

%% Exporting

varname_BLP = ["constant", "coverage", "network","satisfaction","premium", "sigma_coverage","sigma_network","sigma_satisfaction"]';
BLP_estimate = table;
BLP_estimate.name = varname_BLP;
BLP_estimate.coeff = theta_0;
BLP_estimate.se = blpse;
BLP_estimate.tStat = tstat_blp;

save('BLP_estimate.mat','BLP_estimate')

%% Calculate BLP elasticity
ablp = theta1_0(end);
inp.ablp = ablp;

inp.logit = 0;
[elamat_blp] = calela_blp(theta2_0,inp);

%%
elamat_fin_BLP = drvtoela(elamat_blp,inp);
ownpela_blp = diag(elamat_blp);
mcblp = avgp_uniq + (avgshare_uniq./ownpela_blp);
markup_blp_level = avgp_uniq-mcblp;
markup_blp = 100*(avgp_uniq-mcblp)./avgp_uniq; 
inp.mcblp = mcblp;

%% BLP_counterfactual
% Do this counterfactual in terms of the plan level

% To do counterfactual we need function of price which gives the optimal
% price by setting this function equal to zero
% This function is firms' profit maximizing condition w.r.t price
% Given prices, we recalculate 
% (1) fitted share,(2) elasticity using the parameter found 
% and use the backed out marginal costs

% test 1: check that first-order conditions are satisfied in observed prices.
foc_subsidy = @(price) profitFOC(price, theta_0, inp); % function to set to zero
disp('Norm of FOCs in observed prices should be zero. It is:')
foc0 = foc_subsidy(avgp_uniq);
disp(foc0'*foc0)

%% Q1
foc_postsubsidy = @(price) profitFOC(price, theta_0, inp); % function to set to zero

price0 = avgp_uniq; 
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
    'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
    'MaxIterations',25);
[Price_postsubsidy,fval] = fsolve(foc_postsubsidy,price0,options); % starting value = observed prices
%% Q1-1: price chage?
price_pre = avgp_uniq;
perc_ch = 100*(Price_postsubsidy - price_pre)./price_pre;
[~,P_postsubsidy] = foc_postsubsidy(Price_postsubsidy);

%% share change?
sharepost_uniq = zeros(nP,1);

for p =1:nP
    ip = PlanID == p;
    sharepost_uniq(p,:) = mean(P_postsubsidy(ip,:));
end

%% Q1-2: How much will the uninsurance rate change?

share0_bym = zeros(nM,1);
for m = 1:nM
    im = marketID == m;
    share0_bym(m,1) = 1-sum(share(im,1));
end

ogshare_pre = mean(share0_bym);

share0_bym_post = zeros(nM,1);
for m=1:nM
    im = marketID == m;
    share0_bym_post(m,:) = 1-sum(P_postsubsidy(im,:));
end

ogshare_post = mean(share0_bym_post);

%% Q2-Profits

avgshare_uniq_post = zeros(nP,1);
for p = 1:nP
    ip = PlanID ==p;
    avgshare_uniq_post(p,1) = mean(P_postsubsidy(ip,1));
end

profit_pre = (avgp_uniq - mcblp).*avgshare_uniq;
profit_post = (Price_postsubsidy -(mcblp-0.25.*ones(nP,1))).*avgshare_uniq_post;

profit_ch = 100*(profit_post - profit_pre)./profit_pre;

%% Q3-ConsumerSurplus

cs_pre = calcs(avgp_uniq, theta_0,inp);
cs_post = calcs(Price_postsubsidy,theta_0,inp);

cs_ch = -(1/ablp).*(log(cs_post) - log(cs_pre));

%% MC and prodcharac

mcv = mcblp(PlanID);

%%
mc_pcharac = fitlm(prodcharac,mcv);

varname_mc = ["constant", "coverage", "network","satisfaction"]';
mc_estimate = table;
mc_estimate.name = varname_mc;
mc_estimate.coeff = mc_pcharac.Coefficients.Estimate;
mc_estimate.se =mc_pcharac.Coefficients.SE ;
mc_estimate.tStat = mc_pcharac.Coefficients.tStat;
mc_estimate.pValue = mc_pcharac.Coefficients.pValue;


save('mc_estimate.mat','mc_estimate')




