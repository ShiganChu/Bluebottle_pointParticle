function [Lx, Ly, Lz] = cgns_read_part_moment(casename, time)
% CGNS_READ_PART_MOMENT  Read the particle moments from a BLUEBOTTLE-generated
%   CGNS file.
%
%   [Lx, Ly, Lz] = CGNS_READ_PART_MOMENT(CASENAME, TIME) reads the particle
%   body moments from the simulation CASENAME at time TIME. Each of the moment
%   components  is an array representing all of the particles in the simulation
%
%   Example:
%     cgns_read_part_moment('simulation', 3.14159) will read the
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

path = [casename '/output/part-' tt '.cgns'];

lxsol = '/Base/Zone0/Solution/MomentX/ data';
lysol = '/Base/Zone0/Solution/MomentY/ data';
lzsol = '/Base/Zone0/Solution/MomentZ/ data';

Lx = h5read(path, lxsol);
Ly = h5read(path, lysol);
Lz = h5read(path, lzsol);
