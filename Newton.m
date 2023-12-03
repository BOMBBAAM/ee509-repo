function [x_min,f_min,Xk,Fk,Gk,nF,nG,nH,IFLAG] = Newton(x0,epsilon,e_rel,e_abs,itmax)
    % initiate params    
    Xk = [];Fk = []; Gk = [];
    nF = 0;nG = 0;nH = 0;
    IFLAG = -999;
    
    for i = 1:itmax
        [f,gradient,Hessian] = Rosenbrock(x0,3);
        % newton search direction
        s = -inv(Hessian)*gradient;
        % exact line search
        [~, lambda] = golden_2d(x0,x0 + [1.5;1.5],epsilon,itmax);
        
        % update
        x0 = x0 + lambda*s;
        % store
        Xk = [Xk,x0];
        Fk = [Fk;f];
        Gk = [Gk,gradient];
        
        % terminate
        if abs(s .* gradient) <= abs(s .* Gk(:, end)') * e_rel + e_abs
            IFLAG = 1; % converged
            break
        end
    end
    f_min = Fk(end);
    x_min = Xk(:,end);
    
end