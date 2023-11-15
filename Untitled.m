x0=[-1.2,1]';
epsilon = 1e-6;
mu = 1e-4;
eta = 0.25;
itmax = 100;
options = 2;
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(x0,epsilon,mu,eta,itmax,options)

function [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG, nReset] = CG(x0, epsilon, mu, eta, itmax, option)
    g = zeros(2, itmax);
    s = zeros(2, itmax);
    nReset = zeros(itmax, 1);
    xmin = zeros(2, itmax);
    fmin = zeros(itmax, 1);
    Xk = zeros(2, itmax);
    Fk = zeros(itmax, 1);
    Gk = zeros(itmax, 1);
    Lk = zeros(itmax, 1);
    beta = zeros(itmax, 1);
    error = 1;
    delta = 0.1;
    nF = 0;
    nG = 0;
    IFLAG = -99;

    k = 1;
    xmin(:, k) = x0;
    [fmin(k), g(:, k)] = Rosenbrock(x0, option); % pass the 'option' argument
    s(:, k) = -g(:, k);

    while error > epsilon
        Lk(k) = 1;
        xmin(:, k+1) = xmin(:, k) + Lk(k) * s(:, k);
        [fmin(k+1), g(:, k+1)] = Rosenbrock(xmin(:, k+1), option); % pass the 'option' argument

        while (fmin(k+1) > fmin(k) + mu * Lk(k) * (s(:, k)' * g(:, k))) && ...
              (abs(s(:, k)' * g(:, k+1))) > -eta * (s(:, k)' * g(:, k))
            Lk(k) = Lk(k) * delta;
            xmin(:, k+1) = xmin(:, k) + Lk(k) * s(:, k);
            [fmin(k+1), g(:, k+1)] = Rosenbrock(xmin(:, k+1), option);
        end

        if option == 1
            beta(k) = (g(k+1)' * g(k+1)) / (g(k)' * g(k));
        else
            beta(k) = (g(k+1)' * (g(k+1) - g(k))) / (g(k)' * g(k));
        end

        if acos((-s(:, k)' * g(:, k)) / (norm(s(:, k)) * norm(g(:, k)))) > 85 * pi / 180
            beta(k) = 0;
            nReset(k) = 1;
        end

        s(:, k+1) = -g(:, k+1) + beta(k) * s(:, k);
        error = norm(xmin(:, k+1) - xmin(:, k));

        if k > itmax
            disp('not convergence');
            IFLAG = 0;
            return
        end

        if k < itmax
            disp('converged');
            IFLAG = 1;
            for i = 1:k-1
                iter_fmin_xmin_g_IR = [i, fmin(i), xmin(:, i)', g(:, i)', nReset(i)];
                disp(iter_fmin_xmin_g_IR);
                if fmin(i) == 0
                    break
                end
            end
        end

        k = k + 1;
    end
end

function [f, gradient, Hessian] = Rosenbrock(x, options)
    f = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    Hessian = 'NaN';
    gradient = 'NaN';

    if options >= 2
        gradient = [2*x(1) - 400*x(1)*(x(2) - x(1)^2) - 2; -200*x(1)^2 + 200*x(2)];
    end

    if options == 3
        Hessian = [1200*x(1)^2 - 400*x(2) + 2, -400*x(1); -400*x(1), 200];
    end
end

