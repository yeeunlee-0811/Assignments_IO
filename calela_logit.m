function ela_mat_fin = calela_logit(inp)

Krc = inp.Krc;
nP = inp.nP;
nM = inp.nM;
marketID = inp.marketID;
PlanID = inp.PlanID;
al = inp.al;

% use estimated share: P_logit
sigmal = zeros(Krc,1); 
[~,P_logit] = getdelta(sigmal,inp);

ela_mat = zeros(nP,nP,nM);
for i = 1:nP
    for j = 1:nP
        for l = 1:nM
            imi = marketID == l & PlanID == i;
            imj = marketID == l & PlanID == j;
            if any(imi)
                if any(imj)
                    if i == j
                        ela_mat(i,j,l) = al.*(1-P_logit(imi)).*P_logit(imi);
                    else
                        ela_mat(i,j,l) = -al.*P_logit(imj).*P_logit(imi);
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