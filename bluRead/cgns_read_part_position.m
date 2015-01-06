function [x, y, z] = cgns_read_part_position(casename, time)
% CGNS_READ_PART_POSITION Read the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [x, y, z] = CGNS_READ_PART_POSITION(CASENAME, TIME) reads the particle
%   positions from the simulation CASENAME at time TIME. Each of the
%   position components is an array representing all of the particles in
%   the simulation.
%
%   Example:
%     cgns_read_part_position('simulation', 3.14159) will read the
%     appropriate output file located in 'simulation/output

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
 

% path = [casename '/output/part-' tt '.cgns'];
path = [casename '/output/point-' tt '.cgns'];

xsol = '/Base/Zone0/GridCoordinates/CoordinateX/ data';
ysol = '/Base/Zone0/GridCoordinates/CoordinateY/ data';
zsol = '/Base/Zone0/GridCoordinates/CoordinateZ/ data';

x = h5read(path, xsol);
y = h5read(path, ysol);
z = h5read(path, zsol);

