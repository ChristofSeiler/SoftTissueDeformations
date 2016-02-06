function y = piecewise_3_tissues(x)

y = 2.75/100*x.*(x<=100);
y = y + 5.5/100*x.*(x>100 & x<=200);
y = y + 22/100*x.*(x>200);