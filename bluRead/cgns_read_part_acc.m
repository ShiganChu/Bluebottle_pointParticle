function [udot, vdot, wdot] = cgns_read_part_acc(casename, time)
% CGNS_READ_PART_ACC Read the particle accelerations from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [udot, vdot, wdot] = CGNS_READ_PART_ACC(CASENAME, TIME) reads the particle
%   accelerations from the simulation CASENAME at time TIME. Each of the
%   acceleration components is an array representing all of the particles
%   in the simulation.
%
%   Example:
%     cgns_read_part_acc('simulation', 3.14159) will read the appropriate
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

path = [casename '/output/part-' tt '.cgns'];

usol = '/Base/Zone0/Solution/AccelerationX/ data';
vsol = '/Base/Zone0/Solution/AccelerationY/ data';
wsol = '/Base/Zone0/Solution/AccelerationZ/ data';

udot = h5read(path, usol);
vdot = h5read(path, vsol);
wdot = h5read(path, wsol);
