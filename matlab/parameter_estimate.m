clear all;
close all;

% show separat energy
[x y a b0 b1 g] = textread( '../ParameterEstimation/parameterFile.txt', '%f %f %f %f %f %f' );
figure; hold on
xlabel('iterations')
ylabel('energy')
plot(x,y,'--rs')
h = legend('least square difference',1);
set(h,'Interpreter','none')

figure; hold on
xlabel('iterations')
ylabel('parameter')
plot(x,a,'--rs')
h = legend('alpha',1);
set(h,'Interpreter','none')

figure; hold on
xlabel('iterations')
ylabel('parameter')
plot(x,b0,'--rs',x,b1,'--bs',x,g,'--gs')
h = legend('beta0','beta1','gamma',1);
set(h,'Interpreter','none')

yprime = y/max(y);
aprime = a/max(a);
figure; hold on
xlabel('iterations')
ylabel('nomralized values')
plot(x,yprime,'r',x,aprime,'b')
h = legend('least square difference', 'alpha',1);
set(h,'Interpreter','none')

minimumLSD = min(y)