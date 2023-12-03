function [x_min, f_min] = golden_2d(xstart, xstop, epsilon, itmax)
    % golden section params
    tau = 0.5 * (sqrt(5) - 1);
    a1 = xstart(1); b1 = xstart(2);
    a2 = xstop(1); b2 = xstop(2);

    % Initialize points
    x1 = a1 + (1 - tau) * (a2 - a1);
    y1 = b1 + (1 - tau) * (b2 - b1);
    x2 = a1 + tau * (a2 - a1);
    y2 = b1 + tau * (b2 - b1);

    % Define the Rosenbrock function
    func = @(x) Rosenbrock(x,1);

    % Iterate to find the minimum value
    for i = 1:itmax
        f11 = func([x1; y1]); f12 = func([x1; y2]);
        f21 = func([x2; y1]); f22 = func([x2; y2]);

        [fmin, index] = min([f11, f12, f21, f22]);

        if index == 1
            a2 = x2; b2 = y2;
            xmin = [a1;b1];
        elseif index == 2
            a2 = x2; b1 = y1;
            xmin = [a1;b2];
        elseif index == 3
            a1 = x1; b2 = y2;
            xmin = [a2;b1];
        elseif index == 4
            a1 = x1; b1 = y1;
            xmin = [a2;b2];
        end

        % Update points
        x1 = a1 + (1 - tau) * (a2 - a1);
        y1 = b1 + (1 - tau) * (b2 - b1);
        x2 = a1 + tau * (a2 - a1);
        y2 = b1 + tau * (b2 - b1);

        % Check the convergence
        if abs(fmin) < epsilon
            break
        end
    end

    % Update x_min, f_min
    f_min = fmin;
    x_min = xmin;
end
