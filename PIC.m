L = 100;    % domain of solution 0 <= x <= L
N = 20000;  % number of electrons
J = 1000;   % number of grid points
vb = 3;     % beam velocity
dt = 0.1;   % time-step (in inverse plasma frequencies)
tmax = 15;  % simulation run from t = 0 to t = tmax

t = 0;
rng(42);                     % seed the rand # generator
r = L*rand(N,1);             % electron positions
v = double_maxwellian(N,vb); % electron velocities

while (t<=tmax)
    % load r,v into a single vector
    solution_coeffs = [r; v];
    % take a 4th order Runge-Kutta timestep
    k1 = AssembleRHS(solution_coeffs,L,J,N);
    k2 = AssembleRHS(solution_coeffs + 0.5*dt*k1,L,J,N);
    k3 = AssembleRHS(solution_coeffs + 0.5*dt*k2,L,J,N);
    k4 = AssembleRHS(solution_coeffs + dt*k3,L,J,N);
    solution_coeffs = solution_coeffs + dt/6*(k1+2*k2+2*k3+k4);
    % unload solution coefficients
    r = solution_coeffs(1:N);
    v = solution_coeffs(N+1:2*N);
    % make sure all coordinates are in the range 0 to L
    r = r + L*(r<0) - L*(r>L);
    t = t + dt; 
    phase_space(t,dt,r,v);
end
