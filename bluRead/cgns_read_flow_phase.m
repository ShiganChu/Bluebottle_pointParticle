function phase = cgns_read_flow_phase(casename, time)
% CGNS_READ_FLOW_PHASE Read the fluid/solid phase field from a BLUEBOTTLE-
%   generated CGNS file.
%
%   phase = CGNS_READ_FLOW_PHASE(CASENAME, TIME) reads the fluid/solid
%   field from the simulation CASENAME at the time TIME
%
%   Example: 
%     cgns_read_flow_phase('simulation', 3.14159) will read the
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

path = [casename '/output/flow-' tt '.cgns'];

phasesol = '/Base/Zone0/Solution/Phase/ data';

phase = h5read(path, phasesol);
