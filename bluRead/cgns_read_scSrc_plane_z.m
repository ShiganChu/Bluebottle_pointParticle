function p = cgns_read_scalar_plane_z(casename, time, loc)
% CGNS_READ_FLOW_PRES_PLANE_Z Read a plane from the flow velocity field from a
%   BLUEBOTTLE-generated CGNS file.
%
%   p = CGNS_READ_FLOW_PRES_PLANE_Z(CASENAME, TIME, LOC) interpolates the
%   flow pressure field at a plane normal to the z-axis at the location LOC
%   from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_read_flow_pres_plane_z('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     z = 4.2.

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

path = [casename '/output/scalar-' tt '.cgns'];

% get grid for extents
[x y z] = cgns_read_grid(casename);
[sx sy sz] = size(z);
dz = z(1,1,2)-z(1,1,1);
zs = z(1,1,1);
ze = z(1,1,end);
zl = ze - zs;

% check that loc is inside grid
if loc < zs | loc > ze
  disp(sprintf('loc must fall between %f and %f', zs, ze))
  return % quit
end

% get loc index for velocity field interpolation
k = floor((loc - zs)/dz + 0.5);

% get velocity fields
psol = '/Base/Zone0/Solution/SourceScalar/ data';

P = h5read(path, psol);

% do interpolation
if k == 0  % use one-sided extrapolation
  k = 1;
  zz = (k-0.5) * dz + zs;
  p(:,:) = P(:,:,k) + (loc - zz)/dz * (P(:,:,k+1)-P(:,:,k));
elseif k == sz-1  % use one-sided extrapolation
  k = sz-2;
  zz = (k-0.5) * dz + zs;
  p(:,:) = P(:,:,k) + (loc - zz)/dz * (P(:,:,k+1)-P(:,:,k));
else  % use central-difference interpolation
  zz = (k-0.5) * dz + zs;
  p(:,:) = P(:,:,k) + (loc - zz)/dz * (P(:,:,k+1)-P(:,:,k));
end
