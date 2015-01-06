function [sc conv_sc diff_sc scSrc epsp] = cgns_read_scalar(casename, time)
%function [sc conv_sc diff_sc ] = cgns_read_scalar(casename, time)
% CGNS_READ_FLOW_PRES_PLANE_Z Read a plane from the flow velocity field from a
%   BLUEBOTTLE-generated CGNS file.
%
%   p = CGNS_READ_FLOW_PRES_PLANE_Z(CASENAME, TIME, LOC) interpolates the
%   flow pressure field at a plane normal to the z-axis at the location LOC
%   from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_read_flow_pres_plane_z('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     z = 4.2.

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

path = [casename '/output/scalar-' tt '.cgns'];



% get pressure fields
psol = '/Base/Zone0/Solution/Scalar/ data';
convSol = '/Base/Zone0/Solution/ConvScalar/ data';
diffSol = '/Base/Zone0/Solution/DiffScalar/ data';
srcSol = '/Base/Zone0/Solution/SourceScalar/ data';
epspS = '/Base/Zone0/Solution/EPSP/ data';

sc = h5read(path, psol);
conv_sc = h5read(path, convSol);
diff_sc = h5read(path, diffSol);
scSrc = h5read(path, srcSol);
epsp = h5read(path, epspS);


