function [xmin,fmin,IFLAG,nReset] = CG(x0,epsilon,mu,eta,itmax,option)
g = zeros(2,itmax);
s = zeros(2,itmax);
nReset = zeros(itmax);
end