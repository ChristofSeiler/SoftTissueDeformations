clear all;
close all;

direction='vh';
folder='TwoTissueExpansion';

set(0, 'DefaultAxesFontSize', 16);

for i=1:2
    % read data
    [p0 x0 y0] = textread( ['../DVFRegularization/',direction(i),'_plot_0.txt'], '%f %f %f' );
    [p1 x1 y1] = textread( ['../DVFRegularization/',direction(i),'_plot_1.txt'], '%f %f %f' );
    [p2 x2 y2] = textread( ['../DVFRegularization/',direction(i),'_plot_2.txt'], '%f %f %f' );
    [p3 x3 y3] = textread( ['../DVFRegularization/',direction(i),'_plot_3.txt'], '%f %f %f' );
    [p4 x4 y4] = textread( ['../DVFRegularization/',direction(i),'_plot_4.txt'], '%f %f %f' );
    [p5 x5 y5] = textread( ['../DVFRegularization/',direction(i),'_plot_5.txt'], '%f %f %f' );
    [p6 x6 y6] = textread( ['../DVFRegularization/',direction(i),'_plot_6.txt'], '%f %f %f' );
    [p7 x7 y7] = textread( ['../DVFRegularization/',direction(i),'_plot_7.txt'], '%f %f %f' );
    [p8 x8 y8] = textread( ['../DVFRegularization/',direction(i),'_plot_8.txt'], '%f %f %f' );
    [p9 x9 y9] = textread( ['../DVFRegularization/',direction(i),'_plot_9.txt'], '%f %f %f' );

    if(direction(i)=='h') 
        text='horizontal';
        pos=1;
    else
        text='vertical';
        pos=2;
    end
    %xAxis=[text,' profile [pixel]'];
    xAxis='Pixel';
    yAxis='Displacement [pixel]';

    f = figure; hold on
    plot(p0,x0,'-g',p0,y0,'--r', 'LineWidth', 3); title('Level 9', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'9','.eps'],'epsc2');

    f = figure; hold on
    plot(p1,x1,'-g',p1,y1,'--r', 'LineWidth', 3); title('Level 8', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'8','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p2,x2,'-g',p2,y2,'--r', 'LineWidth', 3); title('Level 7', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'7','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p3,x3,'-g',p3,y3,'--r', 'LineWidth', 3); title('Level 6', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'6','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p4,x4,'-g',p4,y4,'--r', 'LineWidth', 3); title('Level 5', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'5','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p5,x5,'-g',p5,y5,'--r', 'LineWidth', 3); title('Level 4', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'4','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p6,x6,'-g',p6,y6,'--r', 'LineWidth', 3); title('Level 3', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'3','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p7,x7,'-g',p7,y7,'--r', 'LineWidth', 3); title('Level 2', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'2','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p8,x8,'-g',p8,y8,'--r', 'LineWidth', 3); title('Level 1', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'1','.eps'],'epsc2');
    
    f = figure; hold on
    plot(p9,x9,'-g',p9,y9,'--r', 'LineWidth', 3); title('Level 0', 'fontweight', 'b'); xlabel(xAxis); ylabel(yAxis); 
    h = legend('x', 'y', pos); set(h,'Interpreter','none');
    saveas(f,['../../../images/',folder,'/',text,'0','.eps'],'epsc2');
end