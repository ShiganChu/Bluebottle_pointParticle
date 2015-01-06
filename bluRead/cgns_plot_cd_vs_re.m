function cgns_plot_cd_vs_re(casename, ts, te, part, nu, rho)
% CGNS_PLOT_CD_VS_RE  Plot the particle drag coefficient vs flow Reynolds number
%   from a BLUEBOTTLE-generated CGNS file.
%
%   CGNS_PLOT_CD_VS_RE(CASENAME, TS, TE, PART) plots the drag coefficient
%   for particle PART vs the flow Reynolds number* from the simulation CASENAME
%   at all times TS <= T <= TE
%
%   Reynolds number is defined using mean flow velocity and particle diameter:
%     Re = U_mean * 2*a / nu
%
%   Example:
%     cgns_plot_cd_vs_re('simulation', 0, 1, 1, 1) will plot the drag
%     coefficient for particle 1 versus flow Reynolds number using NU = 1
%     for the CASENAME = 'simulation' for 0 <= T <= 1

% read part time
[tstr tnum] = cgns_read_flow_time(casename);
n = 1:1:length(tnum);           % create array of indices in tnum

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];          % cut out those indices which do not fall in the time range
ts = tstr(n(1)).time;
te = tstr(n(end)).time;

% calculate the number of timesteps
nt = length(n);

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts);

[np d] = size(x); % np = number of particles

% create data storage locations
CD = zeros(nt, 1);
RE = zeros(nt, 1);

% read flow mean velocity and particle position data
for i = 1:nt
  [u, v, w] = cgns_calc_flow_velocity_mean(casename, tstr(n(i)).time);
  U_mean = (u + v + w) / 3;
  a = cgns_read_part_radius(casename, tstr(n(i)).time);
  RE(i) = U_mean * 2 * a(part) / nu;
  [fx, fy, fz] = cgns_read_part_force(casename, tstr(n(i)).time);
  CD(i) = (fx(part) + fy(part) + fz(part)) / 3;
  CD(i) = CD(i) / (0.5 * rho * U_mean * U_mean * pi * a(part) * a(part));
end

% plot
loglog(RE, CD)
xlabel('Re', 'interpreter', 'latex')
ylabel('$C_D$', 'interpreter', 'latex')
