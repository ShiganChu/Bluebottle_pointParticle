loc = '~/bluebottle/cases/concentration_gradient/';

N = 100;
n = 7;
z = zeros(n*N,1);
ran = [0 1 2 3 4 5 6];

for i = 0:n-1
  loci = sprintf('%s%d', loc, ran(i+1));
  p2 = cgns_read_flow_pres_plane_z(loci, 1.0, 10);
  p1 = cgns_read_flow_pres_plane_z(loci, 1.0, -10);

  dp = mean(mean(p2)) - mean(mean(p1));
  % gradP
  gP(i+1) = dp / 20;

  % histogram of particle z position
  [x y z(i*N+1:i*N+N)] = cgns_read_part_position(loci, 0);
end

hist(z)
