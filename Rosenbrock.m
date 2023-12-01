function [f, gradient, Hessian] = Rosenbrock(x, options)
    f = 100 * (x(2)-x(1)^2)^2 + (1-x(1))^2; % Rosenbrock function
    Hessian = 'NaN'; gradient = 'NaN';
    % calculate gradient
    if options == 3
        Hessian = [1200*x(1)^2 - 400*x(2) + 2,-400*x(1);-400*x(1), 200];
    elseif options >= 2
        gradient = [2*x(1) - 400*x(1)*(x(2) - x(1)^2) - 2;-200*x(1)^2 + 200*x(2)];
    else
        return
    end
end