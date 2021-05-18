% extensive margin
% equations to solve q2 and y in regime 2
function F=eqn_Regime201610(x)
global alpha eta delta s i
global u_y v_q q2 y

y=x;
u_y=y^(-eta);
q2=y*u_y/(1-delta);
v_q=s*q2^(-alpha);

F= u_y*v_q-(1+i);
