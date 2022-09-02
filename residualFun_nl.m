function xi=residualFun_nl(beta,inp)

Xnl = inp.Xnl;
lns_s0 = inp.lns_s0;
jm = inp.jm;

X = [ones(jm,1),Xnl];

xi = lns_s0 - X*beta;
