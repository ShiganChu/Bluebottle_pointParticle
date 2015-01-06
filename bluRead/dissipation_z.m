z = -5.5:1.0:5.5;
Etot = zeros(size(z));

for i = 1:length(z)
  z(i)
  [E k] = cgns_plot_energy_spectrum_z('~/bluebottle/cases/turb_grid/0', 8, 10, z(i));
  Etot(i) = trapz(E);
end

plot(z,Etot,'-')
xlabel('$z$','interpreter','latex')
ylabel('$\langle u_\para^2 \rangle$','interpreter','latex')
