% computing the spreading of influenza
% test for Newton's on the first timestep

% beta transmission, gamma recovery, mu death/birth (replenishment)
beta = 1.00; gamma = 0.20; mu = 0.05;

% initial conditions
N = 500; y02 = 10; y0 = [N-y02 y02 0]';

dt = 1/24; % stepsize for time is h (or dt)
Beta = dt*beta/N; Gamma = dt*gamma; Mu = dt*mu; % for convenience

maxit = 10; tol = 1e-6; % Newton's parameters
fprintf(' k S          I          R          Total     Residual\n');
y = y0;
yinit = y; % initial guess for Newton's
for k = 1:maxit
    % define vector f and its inf norm
    trans_y0 = transpose(y0);
    trans_y = transpose(y);
    f_1 = trans_y0[1] + Mu*N - Mu*trans_y[1] - Beta*trans_y[1]*trans_y[2] - trans_y[1]
    f_2 = trans_y0[2] + Gamma*trans_y[2] - Mu*trans_y[2] + Beta*trans_y[1]*trans_y[2] - trans_y[2]
    f_3 = trans_y0[3] + Gamma*trans_y[2] - Mu*trans_y[3] - trans_y[3]

    f = [f_1, f_2, f_3];
    fnorm = norm(f, inf);
    fprintf('%2d %10.6f %9.6f %9.6f %10.6f %9.2e\n', k-1, y, sum(y), fnorm);
    % stopping criterion
    if fnorm < tol:
        break
    endif
    % define Jacobian matrix
    r_1 = [-Mu -Beta*trans_y[2] - 1, -Beta*trans_y[1], 0];
    r_2 = [Beta*y[2], -Gamma -Mu + Beta * trans_y[1] - 1, 0];
    r_3 = [0, Gamma, -Mu - 1];
    J = [r1;r2;r3];
    J_inv = inv(J);
    % apply Newton's iteration to compute new y
    y = y - J_inv * transpose(f);
end
