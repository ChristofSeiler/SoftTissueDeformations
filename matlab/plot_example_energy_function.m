% plot exponential energy function
x=-5:0.1:5;
y=exp(-x);
f = figure; hold on
title('Minus exponential function')
xlabel('x')
ylabel('y')
plot(x, y, 'LineWidth', 2)
grid on;
saveas(f,'../../../images/eminusx.eps','epsc2');