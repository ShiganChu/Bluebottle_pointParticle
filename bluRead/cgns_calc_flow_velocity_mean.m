function [u, v, w] = cgns_calc_flow_velocity_mean(casename, time)
% CGNS_CALC_FLOW_VELOCITY_MEAN  Read the flow velocity field from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [u, v, w] = CGNS_CALC_FLOW_VELOCITY_MEAN(CASENAME, TIME) calculates the
%   mean flow velocity in each of the three cartesian directions from the
%   simulation CASENAME at time TIME
%   Example:
%     cgns_calc_flow_velocity_mean('simulation', 3.14159) will calculate
%     the mean flow velocity using the appropriate output file located in
%     'simulation/output'

% Sets up the function to take both 'string' and 'double' input for time
if isa(time, 'double') == 1
    tdes = time;
  elseif isa(time, 'char') == 1
    tdes = str2num(time);
end

% Read flow time from file names in casename/output/
[tstr tnum] = cgns_read_flow_time(casename);
n = 1:1:length(tnum);

% find closest time to give time
[c i] = min(abs(tnum - tdes));
time = tstr(i).time;

% read flow velocity
[U, V, W] = cgns_read_flow_vel(casename, time);

% read fluid/solid phase array: phase < 0 ==> fluid cell
phase = cgns_read_flow_phase(casename, time);

% read simulation grid
[X, Y, Z] = cgns_read_grid(casename);

% get grid size (number of cells)
[Xn, Yn, Zn] = size(X);
Xn = Xn - 1;
Yn = Yn - 1;
Zn = Zn - 1;

% get grid extents
Xs = X(1, 1, 1);
Xe = X(Xn,1, 1);
Xl = Xe - Xs;
Ys = Y(1, 1, 1);
Ye = Y(1, Yn,1);
Yl = Ye - Ys;
Zs = Z(1, 1, 1);
Ze = Z(1, 1,Zn);
Zl = Ze - Zs;

% calculate grid spacing
dX = (Xe - Xs) / Xn;
dY = (Ye - Ys) / Yn;
dZ = (Ze - Zs) / Zn;

% get particle radii
a = cgns_read_part_radius(casename, time);
np = size(a);  % number of particles

% calculate volume of fluid phase (subtract volume of particles)
V_part = 0;
for i = 1:np
  V_part = V_part + a(i)*a(i);
end
V_part = V_part * pi;

V_fluid = Xl*Yl*Zl - V_part;

% calculate mean velocity
u = 0;
v = 0;
w = 0;

for k = 1:Zn
  for j = 1:Yn
    for i = 1:Xn
      if phase(i,j,k) < 0
        u = u + U(i,j,k);
        v = v + V(i,j,k);
        w = w + W(i,j,k);
      end
    end
  end
end

% multiply by the size of the cell
% (since all cells are the same size, this factor out of the loop above)
u = u*dX*dY*dZ;
v = v*dX*dY*dZ;
w = w*dX*dY*dZ;

% divide by volume of fluid phase
u = u / V_fluid;
v = v / V_fluid;
w = w / V_fluid;
