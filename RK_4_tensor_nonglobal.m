function [C1] = RK_4_tensor_nonglobal(C0,func,U1,F_tau,t0,t1,tau,A,d)
% This function does a RK-step for a tensor valued ODE with the classical
% Runge Kutta method.

h = t1 - t0;

eval1 = func(C0,F_tau,U1,t0,A,d);
eval2 = func(C0+0.5*h*eval1,F_tau,U1,t0+h/2,A,d);
eval3 = func(C0+0.5*h*eval2,F_tau,U1,t0+h/2,A,d);
eval4 = func(C0+h*eval3,F_tau,U1,t1,A,d);

C1 = C0 + (1/6)*h*(eval1 + 2*eval2 + 2*eval3 + eval4);

end