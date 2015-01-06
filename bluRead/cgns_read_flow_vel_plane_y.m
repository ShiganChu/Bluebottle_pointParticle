function [u, v, w] = cgns_read_flow_vel_plane_y(casename, time, loc)
% CGNS_READ_FLOW_VEL_PLANE_Y Read a plane from the flow velocity field from a
%   BLUEBOTTLE-generated CGNS file.
%
%   [u, v, w] = CGNS_READ_FLOW_VEL_PLANE_Y(CASENAME, TIME, LOC) interpolates the
%   flow velocity field at a plane normal to the y-axis at the location LOC
%   from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_read_flow_vel_plane_y('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     y = 4.2.

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
[sx sy sz] = size(y);
dy = y(1,2,1)-y(1,1,1);
ys = y(1,1,1);
ye = y(1,end,1);
yl = ye - ys;

% check that loc is inside grid
if loc < ys | loc > ye
  disp(sprintf('loc must fall between %f and %f', ys, ye))
  return % quit
end

% get loc index for velocity field interpolation
j = floor((loc - ys)/dy + 0.5);

% get velocity fields
usol = '/Base/Zone0/Solution/VelocityX/ data';
vsol = '/Base/Zone0/Solution/VelocityY/ data';
wsol = '/Base/Zone0/Solution/VelocityZ/ data';

U = h5read(path, usol);
V = h5read(path, vsol);
W = h5read(path, wsol);

% do interpolation
if j == 0  % use one-sided extrapolation
  j = 1;
  yy = (j-0.5) * dy + ys;
  u(:,:) = U(:,j,:) + (loc - yy)/dy * (U(:,j+1,:)-U(:,j,:));
  v(:,:) = V(:,j,:) + (loc - yy)/dy * (V(:,j+1,:)-V(:,j,:));
  w(:,:) = W(:,j,:) + (loc - yy)/dy * (W(:,j+1,:)-W(:,j,:));
elseif j == sy-1  % use one-sided extrapolation
  j = sy-2;
  yy = (j-0.5) * dy + ys;
  u(:,:) = U(:,j,:) + (loc - yy)/dy * (U(:,j+1,:)-U(:,j,:));
  v(:,:) = V(:,j,:) + (loc - yy)/dy * (V(:,j+1,:)-V(:,j,:));
  w(:,:) = W(:,j,:) + (loc - yy)/dy * (W(:,j+1,:)-W(:,j,:));
else  % use central-difference interpolation
  yy = (j-0.5) * dy + ys;
  u(:,:) = U(:,j,:) + (loc - yy)/dy * (U(:,j+1,:)-U(:,j,:));
  v(:,:) = V(:,j,:) + (loc - yy)/dy * (V(:,j+1,:)-V(:,j,:));
  w(:,:) = W(:,j,:) + (loc - yy)/dy * (W(:,j+1,:)-W(:,j,:));
end
