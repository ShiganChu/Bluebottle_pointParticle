function cgns_plot_part_vel(casename, ts, te)
% CGNS_PLOT_PART_VELOCITY  Read the particle velocities from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_VELOCITY(CASENAME, SIGFIGS) plots the particle velocites
%   from the simulation CASENAME between ts-te. Each of the velocity
%   components is an array representing all of the particles in the simulation.
%
%   Example:
%     cgns_plot_part_velocity('simulation',0,1) will plot the velocity history
%     for the CASENAME = 'simulation' between t = 0 and 1 

% read part time
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];
ts = tstr(n(1)).time;
te = tstr(n(end)).time;

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_vel(casename, ts);

nt = length(n); % nt = number of timesteps
[np d] = size(x); % np = number of particles

% create data storage locations
T = zeros(nt, 1);
up = zeros(np, nt);
vp = zeros(np, nt);
wp = zeros(np, nt);

% read particle position data
for i = 1:nt
  [up(:,i), vp(:,i), wp(:,i)] = cgns_read_part_vel(casename, tstr(n(i)).time);
end

% plot
figure
plot(tnum, up)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$u$', 'interpreter', 'latex')

figure
plot(tnum, vp)
xlabel('$t$', 'interpreter','latex')
ylabel('$v$', 'interpreter', 'latex')

figure
plot(tnum, wp)
xlabel('$t$', 'interpreter','latex')
ylabel('$w$', 'interpreter', 'latex')

