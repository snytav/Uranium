function RHS = AssembleRHS( solution_coeffs, L, J )
  r = solution_coeffs(1:N);    v = solution_coeffs(N+1:2*N);
  r = r + L*(r<0) - L*(r>L);
  % Calculate electron number density
  ne = GetDensity( r, L, J );
  % Solve Poisson's equation
  n0 = N/L;
  rho = ne/n0 - 1;
  phi = Poisson1D( rho, L );
  % Calculate electric field
  E = GetElectric( phi, L );
  % equations of motion
  dx = L/J;
  js = floor(r0/dx)+1;
  ys = r0/dx - (js-1);
  js_plus_1 = mod(js,J)+1;
  Efield = E(js).*(1-ys) + E(js_plus_1);
  rdot = v;
  vdot = -Efield;
  RHS = [rdot; vdot];
end