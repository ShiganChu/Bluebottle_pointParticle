function cgns_plot_part_position3(casename, ts, te)
% CGNS_PLOT_PART_POSITION3  Read the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_POSITION3(CASENAME) plots the particle positions from
%   the simulation CASENAME at all times in 3-dimensional space
%   Each of the position components is an array representing all of
%   the particles in the simulation.
%
%   Example:
%     cgns_plot_part_position('simulation', ts, te) will plot the position
%     history for the CASENAME = 'simulation' between ts and te

% find the timestep values
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];
ts = tstr(n(1)).time;
te = tstr(n(end)).time;

% calculate the number of time steps
nt = length(n);

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts);

[np d] = size(x); % np = number of particles

% create data storage locations
X = zeros(np, nt);
Y = zeros(np, nt);
Z = zeros(np, nt);

% read particle position data
for i = 1:nt
  [X(:,i), Y(:,i), Z(:,i)] = cgns_read_part_position(casename, tstr(n(i)).time);
end

% plot
plot3(X,Y,Z)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
zlabel('$z$', 'interpreter', 'latex')
