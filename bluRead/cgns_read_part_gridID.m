function [ii, jj, kk] = cgns_read_part_gridID(casename, time)
% CGNS_READ_PART_VEL  Read the particle velocities from a BLUEBOTTLE-generated
%   CGNS file.
%
%   [u, v, w] = CGNS_READ_PART_VEL(CASENAME, TIME) reads the particle
%   velocities from the simulation CASENAME at time TIME. Each of the
%   velocity components is an array representing all of the particles
%    in the simulation.
%
%   Example:
%     cgns_read_part_vel('simulation', 3.14159) will read the appropriate
%     output file located in 'simulation/output'

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

%path = [casename '/output/part-' tt '.cgns'];
path = [casename '/output/point-' tt '.cgns'];

usol = '/Base/Zone0/Solution/gridIDi/ data';
vsol = '/Base/Zone0/Solution/gridIDj/ data';
wsol = '/Base/Zone0/Solution/gridIDk/ data';

ii = h5read(path, usol);
jj = h5read(path, vsol);
kk = h5read(path, wsol);
