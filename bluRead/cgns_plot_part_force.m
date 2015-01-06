function cgns_plot_part_force(casename, ts, te)
% CGNS_PLOT_PART_FORCE  Plot the particle forces from a BLUEBOTTLE-generated
%   CGNS file.
%
%   CGNS_PLOT_PART_FORCE(CASENAME, TS, TE, DT) plots the particle forces from
%   the simulation CASENAME at all times TS <= T <= TE
% 
%   If TS and TE do not exactly match the flow or part time given in the file
%   name, then the closest time step will be chosen
%
%   Example:
%     cgns_plot_part_force('simulation', 0.0025, 1) will plot the force history
%     for the CASENAME = 'simulation' for 0 <= T <= 1. If the real times 
%     are [0.002 0.0027..... 0.9958, 0.9978, 1.0001] then the force will
%     be plotted for 0.0027 through 0.9978

% read part time
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive 
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];
ts = tstr(n(1)).time;
te = tstr(n(1)).time;

% calculate the number of timesteps
nt = length(n);

% find the number of particles using the initial timestep
[fx, fy, fz] = cgns_read_part_force_total(casename, ts);

[np d] = size(fx); % np = number of particles

% create data storage locations
FX = zeros(np, nt);
FY = zeros(np, nt);
FZ = zeros(np, nt);
FXi = zeros(np, nt);
FYi = zeros(np, nt);
FZi = zeros(np, nt);
FXh = zeros(np, nt);
FYh = zeros(np, nt);
FZh = zeros(np, nt);

% read time and particle force data
for i = 1:nt
  [FX(:,i), FY(:,i), FZ(:,i)] = cgns_read_part_force_total(casename, tstr(n(i)).time);
  [FXi(:,i), FYi(:,i), FZi(:,i)] = cgns_read_part_force_interaction(casename, tstr(n(i)).time);
  [FXh(:,i), FYh(:,i), FZh(:,i)] = cgns_read_part_force_hydro(casename, tstr(n(i)).time);
end

% plot
figure
plot(tnum, FX, tnum, FXi, '--', tnum, FXh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_x$', 'interpreter', 'latex')
figure
plot(tnum, FY, tnum, FYi, '--', tnum, FYh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_y$', 'interpreter', 'latex')
figure
plot(tnum, FZ, tnum, FZi, '--', tnum, FZh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_z$', 'interpreter', 'latex')
