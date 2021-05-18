function ess=eqn_noconnfit1(x)

global year rho vdelta ir

eta=x(1);
s=x(2);
Xstar=x(3);
alpha=x(4);
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

        vz(j)=(1-delta)*qi+yi;
        vGDP(j)=vz(j)+delta*qstar+4*Xstar;
    end
    vzT=vz;
    omega=0.0;
    vtheta=vzT./(vGDP-(1-omega)*4*Xstar);
    vrho=vz./vGDP;
    ess=sum((vrho-rho).^2); % Error Sum of Squares
