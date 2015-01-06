function den = cgns_read_part_radius(casename, time)
% CGNS_READ_PART_RADIUS Read the particle radii from a BLUEBOTTLE-generated
%   CGNS file.
%
%   r = CGNS_READ_PART_RADIUS (CASENAME, TIME) reads the particle radii
%   from the simulation CASENAME at time TIME. The result is an array
%   representing all of the particles in the simulation.
%
%   Example:
%     cgns_read_part_radius('simulation', 3.14159) will read the
%     appropriate output file located in 'simulation/output'

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
path = [casename '/output/point-' tt '.cgns'];
%path = [casename '/output/part-' tt '.cgns'];

rsol = '/Base/Zone0/Solution/Density/ data';

den = h5read(path, rsol);
