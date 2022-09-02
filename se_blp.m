%This calculates variace of blp estimates
function [gmmvar] = se_blp(theta,inp)

jm = inp.jm;
Z = inp.Z_blp;

% m: total number of moments == total number of instruments
m = size(Z,2);

% Gamma
step = 1e-5; % step-length for forward-backward finite differences
K = length(theta);


Gamma = zeros(m,K); 

Dgni = zeros(jm,m,2);
for k=1:K
    for t=1:2 % 1=backward step; 2=forward step;
        theta_eps = theta; % "theta + epsilon" (small perturbation)
        theta_eps(k) = theta_eps(k) + ((-1)^t)*step;
        Dgni(:,:,t) = gi_blp(theta_eps,inp);
    end
    Dgn = sum(Dgni,1);
    Gamma(:,k) = (diff(Dgn,1,3)')/(2*step); % diff between second col. and first col.
end
Gamma = (1/jm).*Gamma;

% V
gni = gi_blp(theta,inp);

gni_squared = zeros(m,m,jm);

for i =1:jm
    gni_squared(:,:,i) = (gni(i,:)')*gni(i,:);
end
V = ((1/jm)^2).*sum(gni_squared,3);


A = (1/jm).*cov(gni);

W = inv(A);




m1 = Gamma'*W*Gamma;
m2 = Gamma'*W*V*W*Gamma;
m3 = inv(m1);
m4 = inv(m2);

% gmmvar0 = inv(m1*(m2\m1));
% gmmvar1 = inv(Gamma'*W*Gamma)*(Gamma'*W*V*W*Gamma)*inv(Gamma'*W*Gamma) ;
gmmvar = inv(Gamma'*(V\Gamma));
end