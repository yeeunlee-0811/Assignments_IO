function gi = gi_blp(theta,inp)

Z = inp.Z_blp;
jm = inp.jm;

m = size(Z,2);

u = residualFun(theta,inp);

gi = zeros(jm,m);

for i =1:m
    gi(:,i)= Z(:,i).*u;
end

end
