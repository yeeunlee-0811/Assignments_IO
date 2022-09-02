function ela_mat_fin = calela_nl(beta,inp)
% calculate dS_k/dp_j 


jm = inp.jm;
nP = inp.nP;
nM = inp.nM;
marketID = inp.marketID;
PlanID = inp.PlanID;
price = inp.premium;
nest = inp.nest;

anl = beta(end,1);
sigma = beta(end-1,1);

% use estimated share: P_nl
P_nl = compshare_nl(beta,price,inp);

Nshare = zeros(jm,1);
for i =1:jm
    in = nest(i,1);
    im = marketID(i,1);
    inm = nest == in & marketID == im; %inm gives index of goods which are in the same nest & in the same market
    Nshare(i,1) = sum(P_nl(inm));
end

% calculate each plan's share in the nest
inNshare = P_nl./Nshare;




ela_mat = zeros(nP,nP,nM);
for i = 1:nP
    for j = 1:nP
        for l = 1:nM
            for n = 1:2
                imi = marketID == l & PlanID == i ;
                imj = marketID == l & PlanID == j;
                if any(imi)
                    if any(imj)
                        if i == j
                            ela_mat(i,j,l) = anl*P_nl(imi)*(1/(1-sigma) - (1/(1-sigma)*inNshare(imi) - P_nl(imi)));
                        else
                            cn1=nest(imi)
                            cn2=nest(imj)
                            if nest(imi) == nest(imj)
                                ela_mat(i,j,l) = -anl.*P_nl(imj)*((sigma/(1-sigma))*inNshare(imi)+P_nl(imi));
                            else
                                ela_mat(i,j,l) = -anl*P_nl(imi).*P_nl(imj);
                            end
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
end
ela_mat_fin = mean(ela_mat,3,'omitnan');



