function [Y1] = RK_4_nonglobal(Y0,tau,func,t0,t1,A,d)
% This is a classical Runge Kutta method of order 4 for time discretization
% of a ODE

h = t1 - t0;

eval1 = func(t0,Y0,A,d);
eval2 = func(t0+h/2,Y0+0.5*h*eval1,A,d);
eval3 = func(t0+h/2,Y0+0.5*h*eval2,A,d);
eval4 = func(t1,Y0+h*eval3,A,d);

Y1 = Y0 + (1/6)*h*(eval1 + 2*eval2 + 2*eval3 + eval4);

end