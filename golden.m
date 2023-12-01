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
    func = @(x) x^4 - 4*sin(x);
    
    % iterate to find minimum value
    for i = 1:itmax
        f1 = func(a); IFunc = IFunc + 1;
        f2 = func(b); IFunc = IFunc + 1;
        if f1 > f2
            a = x1; x1 = x2;
            x2 = b - (1-tau)*(b-a);
            f1 = func(a); IFunc = IFunc + 1;
        else
            b = x2; x2 = x1;
            x1 = a + (1-tau)*(b-a);
            f2 = func(b); IFunc = IFunc + 1;
        end
        
        % check the convergence
        if abs(f1 - f2) < epsilon
            continue
        end
    end
    % update x_min, f_min
    f_min = 0.5*(f1 + f2);
    x_min = (a+b)/2;
    % complete
    IFLAG = 0;
end