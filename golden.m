%% Golden Section
function [x_min,f_min,IFLAG,IFunc] = golden(xstart,epsilon,itmax)
    % IFLAG 0 if success, -999 otherwise
    IFLAG = -999; % initial condition
    % IFunc is he number of function evaluation calls
    IFunc = 0;
    % golden section params
    tau = 0.5*(sqrt(5)-1);
    a = min(xstart); b = max(xstart);
    x1 = a + (1-tau)*(b-a); 
    x2 = b - (1-tau)*(b-a);
    % define function
    func = @sin;
    % iterate to find minimum value
    for i = 1:itmax
        f1 = func(a); f2 = func(b);
        if f1 > f2
            a = x1; x1 = x2;
            x2 = b - (1-tau)*(b-a);
            f1 = func(a);
        else
            b = x2; x2 = x1;
            x1 = a + (1-tau)*(b-a);
            f2 = func(b);
        end
        
        % check the convergence
        if abs(f1 - f2) < epsilon
            f_min = 0.5*(f1 + f2);
            x_min = (a+b)/2;
        end
    end
    % complete
    IFLAG = 0;
end