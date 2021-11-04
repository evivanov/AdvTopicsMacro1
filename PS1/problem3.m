grid = 0:0.1:1;
expon = exp(grid);

intgrid = 0:0.05:1;
interexpon = zeros(length(intgrid), 1);
for i = 1:length(intgrid)
    interexpon(i) = lininterpol(grid, expon, intgrid(i));
end

figure('Name', 'Value Function wo')
plot(intgrid, interexpon);
hold on;
plot(intgrid, exp(intgrid));
legend('interpolated','function','Location','southeast');
saveas(gcf,'ps3','epsc')