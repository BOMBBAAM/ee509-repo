clc;clear;
%% golden
clc;clear;
[xmin, fmin, IFLAG, IFunc] = golden([-1, 2], 1e-6, 100);
disp('Golden Results:');
disp(['xmin: ', num2str(xmin)]);
disp(['fmin: ', num2str(fmin)]);
disp(['IFLAG: ', num2str(IFLAG)]);
disp(['IFunc: ', num2str(IFunc)]);
%% cubic
clc;clear;
[xmin, fmin, IFLAG, IFunc] = cubic([-1, 2], 1e-6, 100);
disp('Cubic Results:')
disp(['xmin: ', num2str(xmin)])
disp(['fmin: ', num2str(fmin)])
disp(['IFLAG: ', num2str(IFLAG)])
disp(['IFunc: ', num2str(IFunc)])
%% newton
clc;clear;
x0 = [0.8;1.5]; epsilon = 1e-6; e_rel = 1e-4; e_abs = 1e-4; itmax = 100;
[x_min,f_min,Xk,Fk,Gk,nF,nG,nH,IFLAG] = Newton(x0,epsilon,e_rel,e_abs,itmax);
disp(x_min)
disp(f_min)