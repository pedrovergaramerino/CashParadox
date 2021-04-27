function ess=eqn_LWfit(x)
% Matching the time series of CIC_GDP ratio in the data with the simulated 
% series in the Lagos-Wright model.

global alpha Xstar year rho s qstar vdelta ir

eta=x(1);
s=x(2);
Xstar=x(3);
alpha=eta;
n=length(year);
vzT=zeros(n,1);
vtheta=vzT;
vz=vzT;
vGDP=vzT;
vrho=vzT;
    for j=1:length(year)
        i=ir(j)/100;
        qi=(s/(1+i))^(1/alpha);
        qstar=s^(1/alpha);
        delta=vdelta(j);
        if delta>1
            vdelta(j)=1;
            delta=1;
        end
        vz(j)=(1-delta)*qi;
        vGDP(j)=(1-delta)*qi+delta*qstar+2*Xstar;
    end
    vzT=vz;
    omega=0.0;
    vtheta=vzT./(vGDP-(1-omega)*2*Xstar);
    vrho=vz./vGDP;
    ess=sum((vrho-rho).^2); % Error Sum of Squares
