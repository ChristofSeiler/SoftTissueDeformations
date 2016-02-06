% plot vector field

clear all;
close all;

%set(0, 'DefaultAxesFontSize', 16);

[x,y] = meshgrid(-4:1:4,-4:1:4);
u=x;
v=y;
f = quiver(x,y,u,v,'LineWidth',2,'Color','k')

title('Vector field', 'fontsize', 16, 'fontweight', 'b');
xlabel('x-Axis', 'fontsize', 16);
ylabel('y-Axis', 'fontsize', 16);

saveas(f,'../../../images/VectorField.eps','epsc2');