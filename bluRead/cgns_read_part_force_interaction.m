function [Fx, Fy, Fz] = cgns_read_part_force_interaction(casename, time)
% CGNS_READ_PART_FORCE_INTERACTION Read the particle forces from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [Fx, Fy, Fz] = CGNS_READ_PART_FORCE_INTERACTION(CASENAME, TIME) reads the
%   particle interaction forces from the simulation CASENAME at time TIME. Each
%   of the force components is an array representing all of the particles in
%   the simulation.
%
%   Example:
%     cgns_read_part_interaction('simulation', 3.14159) will read the
%     appropriate output file located in 'simulation/output'.

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

fxsol = '/Base/Zone0/Solution/InteractionForceX/ data';
fysol = '/Base/Zone0/Solution/InteractionForceY/ data';
fzsol = '/Base/Zone0/Solution/InteractionForceZ/ data';

Fx = h5read(path, fxsol);
Fy = h5read(path, fysol);
Fz = h5read(path, fzsol);
