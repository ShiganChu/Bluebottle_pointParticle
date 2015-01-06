function [x, y, z] = cgns_read_grid(casename)
% CGNS_READ_GRID  Read the discretization grid from a BLUEBOTTLE-generated CGNS
%   file.
%
%   [x, y, z] = CGNS_READ_GRID(CASENAME) reads the flow velocity field from the
%   simulation CASENAME. Each component holds a three-dimensional array
%   whose values provide the value of the corresponding Cartesian vector at
%   each position. For example, the coordinates of a single point in 3D space
%   can be constructed by: x(i,j,k)*i_hat + y(i,j,k)*j_hat + z(i,j,k)*k_hat.
%
%   Example:
%     cgns_read_grid('simulation') will read the appropriate output file located
%     in 'simulation/output'.

path = casename;
path = [path '/output/grid.cgns'];

xsol = '/Base/Zone0/GridCoordinates/CoordinateX/ data';
ysol = '/Base/Zone0/GridCoordinates/CoordinateY/ data';
zsol = '/Base/Zone0/GridCoordinates/CoordinateZ/ data';

x = h5read(path, xsol);
y = h5read(path, ysol);
z = h5read(path, zsol);
