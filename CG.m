function [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(x0,epsilon,mu,eta,itmax,option)
g = zeros(2,itmax);
s = zeros(2,itmax);
nReset = zeros(itmax);
xmin = x0;
fmin = zeros(itmax);
Xk = zeros(2,itmax);
Fk = zeros(itmax);
Lk = zeros(itmax);
nF = 0;
nG = 0;
IFLAG = -99;
[fmin(k), g(:,k)] = Rosenbrock(x0);


end