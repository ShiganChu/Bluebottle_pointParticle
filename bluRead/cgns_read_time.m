function ttime = cgns_read_time(casename, time)
% CGNS_READ_FLOW_VEL  Read the flow velocity field from a BLUEBOTTLE-generated
%   CGNS file.
%
%   [u, v, w] = CGNS_READ_FLOW_VEL(CASENAME, TIME) reads the flow
%   velocity field from the simulation CASENAME at time TIME
%
%   Example:
%     cgns_read_flow_vel('simulation', 3.14159) will read the appropriate
%     output file located in 'simulation/output'

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

tSol = '/Base/Zone0/Etc/Time/ data';
ttime = h5read(path, tSol);


