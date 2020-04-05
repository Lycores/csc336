% computing the spreading of influenza

% beta transmission, gamma recovery, mu death/birth (replenishment)
beta = 1.00; gamma = 0.20; mu = 0.05;

% initial conditions
N = 500; y02 = 10; y0 = [N-y02 y02 0]'; tend = 50;

dt = 1/24; nstep = tend/dt; % stepsize in time and number of time points
Beta = dt*beta/N; Gamma = dt*gamma; Mu = dt*mu;

maxit = 10; tol = 1e-6; % Newton's parameters
y = y0; yi(:, 1) = y;
array = []
for i = 1:nstep % nstep*dt days
    if i >= 5/dt:
        beta = 0.15;
        Beta = dt*beta/N;
    end

    yinit = y; % initial guess for Newton's of the i-th step
    y0 = y;
    array = [array i*dt]
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
        % define Jacobian matrix
        r_1 = [-Mu -Beta*trans_y[2] - 1, -Beta*trans_y[1], 0];
        r_2 = [Beta*y[2], -Gamma -Mu + Beta * trans_y[1] - 1, 0];
        r_3 = [0, Gamma, -Mu - 1];
        J = [r1;r2;r3];
        J_inv = inv(J);
        % apply Newton's iteration to compute new y
        y = y - J_inv * transpose(f);
    end
    
    yi(:, i+1) = y;
    %nit(i) = k;
end

t = (0:nstep)*dt; yn = y;
fprintf('         S      I      R      Total\n');
fprintf('initial: %6.2f %6.2f %6.2f %6.2f\n', y0, sum(y0));
fprintf('end    : %6.2f %6.2f %6.2f %6.2f\n', yn, sum(yn));
fprintf('max infected: max %7.2f\n', max(yi(2, :)));

% do the plot
figure
plot(array, yi[1, :], array, yi[2, :], '--',array, yi[3, :], '-.')
