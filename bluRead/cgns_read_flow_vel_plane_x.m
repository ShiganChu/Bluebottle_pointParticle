function [u, v, w] = cgns_read_flow_vel_plane_x(casename, time, loc)
% CGNS_READ_FLOW_VEL_PLANE_X Read a plane from the flow velocity field from a
%   BLUEBOTTLE-generated CGNS file.
%
%   [u, v, w] = CGNS_READ_FLOW_VEL_PLANE_X(CASENAME, TIME, LOC) interpolates the
%   flow velocity field at a plane normal to the x-axis at the location LOC
%   from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_read_flow_vel_plane_x('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     x = 4.2.

% determine input type of 'time'
if isa(time, 'double') == 1
  % find the directory contents
  path = [casename '/output'];
  od = cd(path);
  contents = dir;
  % take 3rd listing - will be 0.00etc because . and .. exist
  zero = contents(3).name;
  zero = zero(6:end-5);
  sigfigs = length(zero) - 2;
  t_format = sprintf('%%1.%df', sigfigs);
  tt = sprintf(t_format, time);
  cd(od);
elseif isa(time, 'char') == 1
  tt = time;
end

path = [casename '/output/flow-' tt '.cgns'];

% get grid for extents
[x y z] = cgns_read_grid(casename);
[sx sy sz] = size(x);
dx = x(2,1,1)-x(1,1,1);
xs = x(1,1,1);
xe = x(end,1,1);
xl = xe - xs;

% check that loc is inside grid
if loc < xs | loc > xe
  disp(sprintf('loc must fall between %f and %f', xs, xe))
  return % quit
end

% get loc index for velocity field interpolation
i = floor((loc - xs)/dx + 0.5);

% get velocity fields
usol = '/Base/Zone0/Solution/VelocityX/ data';
vsol = '/Base/Zone0/Solution/VelocityY/ data';
wsol = '/Base/Zone0/Solution/VelocityZ/ data';

U = h5read(path, usol);
V = h5read(path, vsol);
W = h5read(path, wsol);

% do interpolation
if i == 0  % use one-sided extrapolation
  i = 1;
  xx = (i-0.5) * dx + xs;
  u(:,:) = U(i,:,:) + (loc - xx)/dx * (U(i+1,:,:)-U(i,:,:));
  v(:,:) = V(i,:,:) + (loc - xx)/dx * (V(i+1,:,:)-V(i,:,:));
  w(:,:) = W(i,:,:) + (loc - xx)/dx * (W(i+1,:,:)-W(i,:,:));
elseif i == sx-1  % use one-sided extrapolation
  i = sx-2;
  xx = (i-0.5) * dx + xs;
  u(:,:) = U(i,:,:) + (loc - xx)/dx * (U(i+1,:,:)-U(i,:,:));
  v(:,:) = V(i,:,:) + (loc - xx)/dx * (V(i+1,:,:)-V(i,:,:));
  w(:,:) = W(i,:,:) + (loc - xx)/dx * (W(i+1,:,:)-W(i,:,:));
else  % use central-difference interpolation
  xx = (i-0.5) * dx + xs;
  u(:,:) = U(i,:,:) + (loc - xx)/dx * (U(i+1,:,:)-U(i,:,:));
  v(:,:) = V(i,:,:) + (loc - xx)/dx * (V(i+1,:,:)-V(i,:,:));
  w(:,:) = W(i,:,:) + (loc - xx)/dx * (W(i+1,:,:)-W(i,:,:));
end
