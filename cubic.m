function [x_min,f_min,IFLAG,IFunc] = cubic(xstart,epsilon,itmax)
    % IFLAG 0 if success, -999 otherwise
    IFLAG = -999; % initial condition
    % IFunc is he number of function evaluation calls
    IFunc = 0;
    % Note: f'(A) < 0 and f'(B) > 0 and A < B
    % initialize cubic params
    A = min(xstart); B = max(xstart);
    % define function
    f = @(x) x^4 - 4*sin(x); IFunc = IFunc + 1;
    df = @(x) 4*x^3 - 4*cos(x);
    
    for i = 1:itmax
        Z = 3*((f(A) - f(B))/(B-A)) + df(A) + df(B); IFunc = IFunc + 2;
        Q = sqrt(Z^2 - df(A)*df(B));
        lambda = A + (df(A) + Z + Q)*(B-A)/(df(A) + df(B) + 2*Z); 

        if f(A) > f(lambda)
            A = lambda;
        else
            B = lambda;
        end
        
        % termination
        if abs(df(lambda)) <= epsilon || ...
                (abs(A - lambda) <= epsilon && abs(B - lambda) <= epsilon)
            break
        end
    end
    % update x_min,f_min
    x_min = lambda;
    f_min = f(lambda); IFunc = IFunc + 1;
    IFLAG = 0;
end