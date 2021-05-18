function ess=eqn_noconnfit2(x)
% NCF model with unregistered activities

global year rho vdelta ir

eta=x(1);
s=x(2);
Xstar=x(3);
alpha=eta;
% alpha=x(4);
n=length(year);
vzT=zeros(n,1);
vtheta=vzT;
vz=vzT;
vGDP=vzT;
vrho=vzT;
    for j=1:length(year)
        i=ir(j)/100;
        qi=(s/(1+i))^(1/alpha);
        yi=(1+i)^(-1/eta);
        qstar=s^(1/alpha);
        if vdelta(j)>1
            vdelta(j)=1;
        elseif vdelta(j)<0
            vdelta(j)=0;
        end
        delta=vdelta(j);
        vzT(j)=(1-delta)*qi;
        vz(j)=vzT(j)+yi;
        vGDP(j)=vzT(j)+delta*qstar+4*Xstar;
    end
    vtheta=vzT./(vGDP-4*Xstar);
    vrho=vz./vGDP;
    ess=sum((vrho-rho).^2); % Error Sum of Squares
