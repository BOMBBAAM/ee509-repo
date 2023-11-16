clc; clear;

% Define parameters
x0 = [-1.2, 1]';
epsilon = 1e-6;
mu = 1e-4;
eta = 0.25;
itmax = 100;
options_FR = 1; % Fletcher-Reeves
options_PR = 2; % Polak-Ribiére

% Run Fletcher-Reeves method
disp('Running Fletcher-Reeves Method:');
[xmin_FR, fmin_FR, ~, ~, ~, ~, nF_FR, nG_FR, IFLAG_FR, nReset_FR] = CG(x0, epsilon, mu, eta, itmax, options_FR);

% Run Polak-Ribiére method
disp('Running Polak-Ribiére Method:');
[xmin_PR, fmin_PR, ~, ~, ~, ~, nF_PR, nG_PR, IFLAG_PR, nReset_PR] = CG(x0, epsilon, mu, eta, itmax, options_PR);

function [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG, nReset] = CG(x0, epsilon, mu, eta, itmax, option)
    % Initialize variables
    g = zeros(length(x0), itmax);
    s = zeros(length(x0), itmax);
    nReset = zeros(itmax, 1);
    xmin = zeros(length(x0), itmax);
    fmin = zeros(itmax, 1);
    Xk = zeros(length(x0), itmax);
    Fk = zeros(itmax, 1);
    Gk = zeros(length(x0), itmax);
    Lk = zeros(itmax, 1);
    beta = zeros(itmax, 1);
    error = 1;
    delta = 0.1;
    nF = 0;
    nG = 0;
    IFLAG = -99;

    k = 1;
    xmin(:, k) = x0;
    [fmin(k), g(:, k)] = Rosenbrock(x0, 2); % Initial function value and gradient
    s(:, k) = -g(:, k);

    % Main optimization loop
    while error > epsilon
        Lk(k) = 1;
        xmin(:, k+1) = xmin(:, k) + Lk(k) * s(:, k);
        [fmin(k+1), g(:, k+1)] = Rosenbrock(xmin(:, k+1), 2);
        nF = nF + 1; % Increment the function evaluation count
        nG = nG + 1; % Increment the gradient evaluation count

        % Line search using Armijo's and Strong Wolfe's conditions
        while (fmin(k+1) > fmin(k) + mu * Lk(k) * (s(:, k)' * g(:, k))) && ...
              (abs(s(:, k)' * g(:, k+1))) > -eta * (s(:, k)' * g(:, k))
            Lk(k) = Lk(k) * delta;
            nF = nF + 1; % Increment the function evaluation count
            xmin(:, k+1) = xmin(:, k) + Lk(k) * s(:, k);
            [fmin(k+1), g(:, k+1)] = Rosenbrock(xmin(:, k+1), 2);
        end


        % Update direction using specified beta formula
        if option == 1 % Fletcher-Reeves
            beta(k) = (g(:, k+1)' * g(:, k+1)) / (g(:, k)' * g(:, k));
        else           % Polak-Ribiere
            beta(k) = (g(:, k+1)' * (g(:, k+1) - g(:, k))) / (g(:, k)' * g(:, k));
        end

        % Watch-dog scheme to reset the search direction
        angle_condition = acos((-s(:, k)' * g(:, k)) / (norm(s(:, k)) * norm(g(:, k))));
        if angle_condition > 85 * pi / 180 || (s(:, k)' * g(:, k+1)) >= 0
            beta(k) = 0;
            nReset(k) = 1;
        end

        % Update search direction
        s(:, k+1) = -g(:, k+1) + beta(k) * s(:, k);
        error = norm(xmin(:, k+1) - xmin(:, k));
        
        % Display information about each iteration
        disp(['Iteration ', num2str(k), ':']);
        disp(['   xk: ', num2str(xmin(:, k+1)')]);
        disp(['   f(xk): ', num2str(fmin(k+1))]);
        disp(['   ∇f(xk): ', num2str(g(:, k+1)')]);
        disp(['   λk: ', num2str(Lk(k))]);
        disp(['   nF: ', num2str(nF)]);
        disp(['   nG: ', num2str(nG)]);
        disp(['   Restart: ', num2str(nReset(k))]);
        disp('------------------------------------');

        % Check for maximum iterations
        if k >= itmax
            disp('Maximum iterations reached. Optimization did not converge.');
            IFLAG = 0;
            return
        end

        k = k + 1;
    end

    IFLAG = 1;
end

function [f, gradient] = Rosenbrock(x, options)
    % Rosenbrock function
    f = 100 * sum((x(2:end) - x(1:end-1).^2).^2) + sum((1 - x(1:end-1)).^2);
    gradient = 'NaN';

    % Calculate gradient if requested
    if options >= 2
        gradient = zeros(length(x), 1);
        gradient(1) = -400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1));
        gradient(end) = 200 * (x(end) - x(end-1)^2);
        gradient(2:end-1) = 200 * (x(2:end-1) - x(1:end-2).^2) - 400 * x(2:end-1) .* (x(3:end) - x(2:end-1).^2);
    end
end

