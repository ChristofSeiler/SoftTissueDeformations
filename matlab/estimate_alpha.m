clear all;
close all;

% show separat energy
[x t p o] = textread( '../DVFRegularization/energyFile.txt', '%f %f %f %f' );
figure; hold on
xlabel('iterations')
ylabel('energy')
plot(x,t,'--rs',x,p,'--bs',x,o,'--gs')
h = legend('combined energy','prior energy','observation energy',1);
set(h,'Interpreter','none')

% diff of energies
[x t1 p1 o1] = textread( '../DVFRegularization/energyFile_1equ.txt', '%f %f %f %f' );
[x t2 p2 o2] = textread( '../DVFRegularization/energyFile_3equ.txt', '%f %f %f %f' );
figure; hold on
xlabel('iterations')
ylabel('energy')
plot(x,t2-t1,'--rs',x,p2-p1,'--bs',x,o2-o1,'--gs')
h = legend('combined energy difference','prior energy difference','observation energy difference',2);
set(h,'Interpreter','none')
