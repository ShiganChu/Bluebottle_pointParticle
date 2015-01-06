function cgns_plot_pmom_vs_gradp(casename, ts, te, dt, part, nu, rho, P, dP)
% CGNS_PLOT_CD_VS_RE  Plot the particle drag coefficient vs flow Reynolds number
%   from a BLUEBOTTLE-generated CGNS file.
%
%   CGNS_PLOT_CD_VS_RE(CASENAME, TS, TE, DT, PART) plots the drag coefficient
%   for particle PART vs the flow Reynolds number* from the simulation CASENAME
%   at all times TS <= T <= TE, with DT matching the output time interval set in
%   record.config.
%
%   Reynolds number is defined using mean flow velocity and particle diameter:
%     Re = U_mean * 2*a / nu
%
%   Example:
%     cgns_plot_cd_vs_re('simulation', 0, 1, 0.01, 1, 1) will plot the drag
%     coefficient for particle 1 versus flow Reynolds number using NU = 1
%     for the CASENAME = 'simulation' for 0 <= T <= 1 outputted at an interval
%     of 0.01.

% calculate the number of timesteps
nt = ceil((te - ts) / dt + 1); % nt = number of timesteps

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts, dt);

[np d] = size(x); % np = number of particles

% create data storage locations
T = zeros(nt, 1);
pmom = zeros(nt, 1);
gradp = zeros(nt, 1);

% read flow mean velocity and particle position data
for i = 1:nt
  time = ts + (i-1)*dt;
  T(i) = cgns_read_part_time(casename, time, dt);
  gradp(i) = abs(dP*T(i));
  if gradp(i) > abs(P)
    gradp(i) = abs(P);
  end
  [fx, fy, fz] = cgns_read_part_force(casename, time, dt);
  pmom(i) = fz(part);%sqrt(fx(part)*fx(part) + fy(part)*fy(part) + fz(part)*fz(part));
  pmom(i) = 1 - pmom(i) / 64 / gradp(i)/(1 - 4/3*pi/64);
end

pmom(end)

% plot
plot(T, pmom)
xlabel('Re', 'interpreter', 'latex')
ylabel('$C_D$', 'interpreter', 'latex')
