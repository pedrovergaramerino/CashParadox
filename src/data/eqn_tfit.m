function ess=eqn_tfit(x)
% Matching the time series of CIC_GDP ratio in the data with the simulated 
% series in the full model

global alpha eta Xstar year ir rho
global s ystar qstar vdelta
global i delta u_y v_q q2 y

eta=x(1);
s=x(2);
Xstar=x(3);
alpha=eta;
% alpha=x(4);
n=length(year);
vq2=zeros(n,1);
vy=vq2;
vz2=vq2;
vz=vz2;
vzB=vq2;
vregime=vq2;
vGDP=vq2;
options=optimset('MaxFunEvals',50000 ,'Display', 'off', 'MaxIter',5000);   % Option to display output

for j=1:length(year)
    i=ir(j)/100;
    qi=(s/(1+i))^(1/alpha);
    qstar=s^(1/alpha);
    yi=(1+i)^(-1/eta);
    delta=vdelta(j);
    if delta>1
        vdelta(j)=1;
        delta=1;
    end
    delta12=1-ystar/qi;     % Cutoff of Regimes 1 and 2
    delta23=1-yi*(1+i)/qstar;   % Cutoff of Regimes 2 and 3
    
    if delta12<0
        NoRegime1=1;
    end
    if delta<=delta12       % Regime 1
        vq2(j)=qi;
        vy(j)=ystar;
        vz2(j)=qi;
        vzB(j)=0;
        vregime(j)=1;
    elseif delta>=delta23   % Regime 3
        vq2(j)=qstar;
        vy(j)=yi;
        vz2(j)=qstar/(1+i);
        vzB(j)=yi-(1-delta)*qstar/(1+i);
        vregime(j)=3;
    else                    % Regime 2
        x0=ystar;
        [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
        vq2(j)=q2;
        vy(j)=y;
        vz2(j)=y/(1-delta);
        vzB(j)=0;
        vregime(j)=2;
    end
end


vzA=(1-vdelta).*vz2;
vz=vzA+vzB;
vzT=(1-vdelta).*vz2+vy;
vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+vy+3*Xstar; % count CM
vrho=vz./vGDP;
ess=sum((vrho-rho).^2); % Error Sum of Squares