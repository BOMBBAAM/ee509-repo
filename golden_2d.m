function [x_min,f_min,IFLAG,IFunc] = golden_2d(xstart,xstop,epsilon,itmax)
    % IFLAG 0 if success, -999 otherwise
    IFLAG = -999; % initial condition
    % IFunc is he number of function evaluation calls
    IFunc = 0;
    
    % golden section params
    tau = 0.5*(sqrt(5)-1);
    a1 = xstart(1); b1 = xstart(2);
    a2 = xstop(1); b2 = xstop(2);
    % (x1,y1)
    x1 = a1 + (1-tau)*(a2-a1);
    y1 = b1 + (1-tau)*(b2-b1);
    % (x2,y2)
    x2 = a2 - (1-tau)*(a2-a1);
    y2 = b2 - (1-tau)*(b2-b1);
    
    % define function
    func = @(x,option) Rosenbrock(x,option);
    
    % iterate to find minimum value
    for i = 1:itmax
        f11 = func([x1,y1]); f12 = func([x1,y2]); 
        f21 = func([x2,y1]); f22 = func([x2,y2]);
        fmin = min([f11;f12;f21;f22]);
        % update range
        if fmin = f11
            a2 = x2; b2 = y2;
        elseif fmin = f12
            a2 = x2; b1 = y1;
        elseif fmin = f21
            a1 = x1; b2 = y2;
        elseif fmin = f22
            a1 = x1; b1 = y2;
        end
        
        % update (x1,y1)
        x1 = a1 + (1-tau)*(a2-a1);
        y1 = b1 + (1-tau)*(b2-b1);
        % update (x2,y2)
        x2 = a2 - (1-tau)*(a2-a1);
        y2 = b2 - (1-tau)*(b2-b1);
        
        % check the convergence
        if sqrt((a2-a1)^2 + (b2-b1)^2) < epsilon
            break
        end
    end
    % update x_min, f_min
    f_min = 0.5*(f1 + f2);
    x_min = (a+b)/2;
    % complete
    IFLAG = 0;
end