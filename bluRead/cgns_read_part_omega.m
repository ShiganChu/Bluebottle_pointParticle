function [Ox, Oy, Oz] = cgns_read_part_omega(casename, time)
% CGNS_READ_PART_OMEGA  Read the particle angular velocities from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [Ox, Oy, Oz] = CGNS_READ_PART_OMEGA(CASENAME, TIME) reads the particle
%   angular velocities from the simulation CASENAME at time TIME. Each of
%   the angular  velocity components is an array representing all of the
%   particles in the  simulation.
%
%   Example:
%     cgns_read_part_omega('simulation', 3.14159) will read the
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

oxsol = '/Base/Zone0/Solution/AngularVelocityX/ data';
oysol = '/Base/Zone0/Solution/AngularVelocityY/ data';
ozsol = '/Base/Zone0/Solution/AngularVelocityZ/ data';

Ox = h5read(path, oxsol);
Oy = h5read(path, oysol);
Oz = h5read(path, ozsol);

