function cgns_plot_part_position(casename, ts, te)
% CGNS_PLOT_PART_POSITION  Plot the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_POSITION(CASENAME, TS, TE) plots the particle positions
%   from the simulation CASENAME at all times TS <= T <= TE
%
%   Example:
%     cgns_plot_part_position('simulation', 0, 1) will plot the position
%     history for the CASENAME = 'simulation' for 0 <= T <= 1

% read part time
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];
ts = tstr(n(1)).time;
te = tstr(n(end)).time;


% calculate the number of timesteps
nt = length(n);

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts);

[np d] = size(x); % np = number of particles

% create data storage locations
X = zeros(np, nt);
Y = zeros(np, nt);
Z = zeros(np, nt);

% read time and particle position data
for i = 1:nt
  [X(:,i), Y(:,i), Z(:,i)] = cgns_read_part_position(casename, tstr(n(i)).time);
end

% plot
figure
plot(tnum, X)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')

figure
plot(tnum, Y)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')

figure
plot(tnum, Z)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
