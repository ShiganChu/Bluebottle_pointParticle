function cgns_plot_part_position(casename, ts, te, dt)
% CGNS_PLOT_PART_POSITION  Plot the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_POSITION(CASENAME, TS, TE, DT) plots the particle positions
%   from the simulation CASENAME at all times TS <= T <= TE, with DT matching
%   the output time interval set in record.config.
%
%   Example:
%     cgns_plot_part_position('simulation', 0, 1, 0.01) will plot the position
%     history for the CASENAME = 'simulation' for 0 <= T <= 1 outputted at an
%     interval of 0.01.

% calculate the number of timesteps
nt = ceil((te - ts) / dt + 1); % nt = number of timesteps

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts, dt);
[u, v, w] = cgns_read_part_velocity(casename, ts, dt);
[fx, fy, fz] = cgns_read_part_force(casename, ts, dt);
[ox, oy, oz] = cgns_read_part_omega(casename, ts, dt);
[lx, ly, lz] = cgns_read_part_moment(casename, ts, dt);

[np d] = size(x); % np = number of particles

% create data storage locations
T = zeros(nt, 1);
X = zeros(nt, np);
Y = zeros(nt, np);
Z = zeros(nt, np);
U = zeros(nt, np);
V = zeros(nt, np);
W = zeros(nt, np);
Fx = zeros(nt, np);
Fy = zeros(nt, np);
Fz = zeros(nt, np);
Ox = zeros(nt, np);
Oy = zeros(nt, np);
Oz = zeros(nt, np);
Lx = zeros(nt, np);
Ly = zeros(nt, np);
Lz = zeros(nt, np);

% read time and particle position data
for i = 1:nt
  time = ts + (i-1)*dt;
  T(i) = cgns_read_part_time(casename, time, dt);
  [X(i,:), Y(i,:), Z(i,:)] = cgns_read_part_position(casename, time, dt);
  [U(i,:), V(i,:), W(i,:)] = cgns_read_part_velocity(casename, time, dt);
  [Fx(i,:), Fy(i,:), Fz(i,:)] = cgns_read_part_force(casename, time, dt);
  [Ox(i,:), Oy(i,:), Oz(i,:)] = cgns_read_part_omega(casename, time, dt);
  [Lx(i,:), Ly(i,:), Lz(i,:)] = cgns_read_part_moment(casename, time, dt);
end

% plot
subplot(5,3,1)
plot(T, X)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')
subplot(5,3,2)
plot(T, Y)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
subplot(5,3,3)
plot(T, Z)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
subplot(5,3,4)
plot(T, U)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$u$', 'interpreter', 'latex')
subplot(5,3,5)
plot(T, V)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$v$', 'interpreter', 'latex')
subplot(5,3,6)
plot(T, W)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$w$', 'interpreter', 'latex')
subplot(5,3,7)
plot(T, Fx)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_x$', 'interpreter', 'latex')
subplot(5,3,8)
plot(T, Fy)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_y$', 'interpreter', 'latex')
subplot(5,3,9)
plot(T, Fz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_z$', 'interpreter', 'latex')
subplot(5,3,10)
plot(T, Ox)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_x$', 'interpreter', 'latex')
subplot(5,3,11)
plot(T, Oy)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_y$', 'interpreter', 'latex')
subplot(5,3,12)
plot(T, Oz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_z$', 'interpreter', 'latex')
subplot(5,3,13)
plot(T, Lx)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_x$', 'interpreter', 'latex')
subplot(5,3,14)
plot(T, Ly)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_y$', 'interpreter', 'latex')
subplot(5,3,15)
plot(T, Lz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_z$', 'interpreter', 'latex')
