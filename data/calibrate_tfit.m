% extensive margin
% count CM - all consume Xstar

% v(q)=s*(q^(1-alpha)/(1-alpha));
% u(y)=y^(1-eta)/(1-eta);
% count CM

clear;

global alpha eta i rhoTarget Xstar delta omega shareZBTarget
global s ystar yi qstar qi year ir i1
global rho u_y v_q q2 y vdelta
global NoRegime1 flag
fig=1;  % Control on generating different results;
        % "=0", shows each model for each country under each interest rate, i.e.,
        %       displays each result one by one.
        % "=1", displays Figures 5(a) to 5(d), A.2(a) to A.2(d) in the paper
        %       by comparing three models: (1) standard model, (2) NCF
        %       model, (3) full model.
        % "=2", generates Figure D.2 in the paper by comparing different
        %       interest rates.

irflag=1;   % Selection of inflation expectations
            % "=1" (default), corporate bond rate
            % "=2", DLM of actual inflation
            % "=3", Term structure based

flag=3; % Country selection;
        % 0: Canada
        % 1: US
        % 2: Australia
        % 3: UK
model=1;    % 0: Lagos-Wright Model (Standard model)
            % 1: JS model (Full model)
            % 2: JS model with different eta and alpha
            % 3: Model shut down type B agents, alpha/=eta
            % 4: Model shut down type B agents, alpha=eta
            % 5: Underground economy JS model (Full-UA)
            % 6: Underground economy NCF model (NCF-UA)

ystar=1;
set(gcf,'DefaultLineLineWidth',1.0);
if (flag==2)
    sn=20;
elseif (flag==3)
    sn=23;
else
    sn=30;
end

if fig==0
    
    if irflag==1
        load data;
        if flag==0 % Canada
            year=can(:,1);
            ir=can(:,2);
            ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
            yt=[1960,1980,1999,2005,2012]';
            rho=can(:,3);
        elseif flag==1 % US
            year=usa(:,1);
            ir=usa(:,2);
            ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
            yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
            rho=usa(:,3);
        elseif flag==2 % Australia
            year=aus(:,1);
            ir=aus(:,2);
            ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                0.90064194]';
            yt=[1965,1994:2011]';
            rho=aus(:,3);
        else  % UK
            uk(1:4,:)=[];
            year=uk(:,1);
            ir=uk(:,2);
            ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
            yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
            rho=uk(:,3);
        end
    else
        load ir
        if flag==0 % Canada
            yr0=ir_can(:,1);
            if irflag==2
                ir0=ir_can(:,3);    % DLM of actual inflation
            else
                ir0=ir_can(:,4);     % Term structure based
            end
            rho0=ir_can(:,6);
            idx=(isnan(ir0))|(isnan(rho0));
            year=yr0(~idx);
            ir=ir0(~idx)*100;
            rho=rho0(~idx);
            ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
            yt=[1960,1980,1999,2005,2012]';
        elseif flag==1 % US
            yr0=ir_usa(:,1);
            if irflag==2
                ir0=ir_usa(:,3);    % DLM of actual inflation
            else
                ir0=ir_usa(:,4);    % Term structure based
            end
            rho0=ir_usa(:,7);
            idx=(isnan(ir0))|(isnan(rho0));
            year=yr0(~idx);
            ir=ir0(~idx)*100;
            rho=rho0(~idx);
            ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
            yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
        elseif flag==2 % Australia
            yr0=ir_aus(:,1);
            if irflag==2
                ir0=ir_aus(:,3);    % DLM of actual inflation
            else
                ir0=ir_aus(:,4);    % Term structure based
            end
            rho0=ir_aus(:,7);
            idx=(isnan(ir0))|(isnan(rho0));
            year=yr0(~idx);
            ir=ir0(~idx)*100;
            rho=rho0(~idx);
            ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                0.90064194]';
            yt=[1960,1994:2011]';
        else  % UK
            yr0=ir_uk(:,1);
            if irflag==2
                ir0=ir_uk(:,3);    % DLM of actual inflation
            else
                ir0=ir_uk(:,4);    % Term structure based
            end
            rho0=ir_uk(:,6);
            idx=(isnan(ir0))|(isnan(rho0));
            year=yr0(~idx);
            ir=ir0(~idx)*100;
            rho=rho0(~idx);
            year(1:4)=[];
            ir(1:4)=[];
            rho(1:4)=[];
            ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
            yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
        end
    end
    % Use the data to find the fitted curve of delta
    [dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
    vdelta=feval(dt,year); % evaluate delta at each year between 1960 to 2014
    n=length(year);
    vq2=zeros(n,1);
    vy=vq2;
    vz2=vq2;
    vz=vz2;
    vzB=vq2;
    vregime=vq2;
    vGDP=vq2;
    vrho=vq2;
    vtheta=vq2;
    vzT=vq2;
    
    options=optimset('MaxFunEvals',50000 ,'Display', 'off', 'MaxIter',5000);
    nonlcon=[];
    fmcopt=optimoptions(@fmincon,'MaxFunEvals',50000,'MaxIter',5000);% Option to display output
    if model==1
        if flag==3
            x0=[0.9 2 10]; % [eta s Xstar] initial guess
        else
            x0=[0.9,1.5,5];
        end
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0,1.01,0];
        ub=[0.999999,inf,inf];
        
        [x,ess,exitflag] = fmincon(@eqn_tfit,x0,A,b,Aeq,beq,lb,ub);
        eta=x(1); s=x(2); Xstar=x(3); alpha=eta;
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
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
            delta12=1-ystar/qi; % Cutoff of Regimes 1 and 2
            delta23=1-yi*(1+i)/qstar; % Cutoff of Regimes 2 and 3
            
            if delta12<0
                NoRegime1=1;
            end
            if delta<=delta12   % Regime 1
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
            else                % Regime 2
                x0=ystar;
                [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
                vq2(j)=q2;
                vy(j)=y;
                vz2(j)=y/(1-delta);
                vzB(j)=0;
                x0=ystar;
                vregime(j)=2;
            end
        end
        
        
        vzA=(1-vdelta).*vz2;
        vz=vzA+vzB;
        vzT=(1-vdelta).*vz2+vy;
        vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+vy+3*Xstar; % count CM
        
        vvelocity=vzT./vz;
        vtheta=vzT./(vGDP-3*Xstar);
        vrho=vz./vGDP;
        
    elseif model==0
        if flag==2
            x0=[0.5,1.5,15];
        else
            x0=[0.9,1.5,10.0];
        end
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0,1,0];
        ub=[0.99999,1,inf];
        [x,ess,exitflag] = fmincon(@eqn_LWfit,x0,A,b,Aeq,beq,lb,ub);
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
        
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
        
    elseif model==2
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0.0,1.0015,0,0.0];
        ub=[0.99999,inf,inf,0.99999];
        if flag==0
            x0=[0.5,1.5,15,0.2];
        elseif flag==3
            x0=[0.956,1.89,8.857,0.956];
            A=[]; b=[]; Aeq=[]; beq=[];
            lb=[0.0,1,0,0.9];
            ub=[0.99999,inf,inf,0.99999];
            
        else
            x0=[0.5,2,10,0.9];
        end
        
        [x,ess,exitflag] = fmincon(@eqn_tfit,x0,A,b,Aeq,beq,lb,ub);
        eta=x(1); s=x(2); Xstar=x(3); alpha=x(4);
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
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
            delta12=1-ystar/qi;
            delta23=1-yi*(1+i)/qstar;
            
            if delta12<0
                NoRegime1=1;
            end
            if delta<=delta12
                vq2(j)=qi;
                vy(j)=ystar;
                vz2(j)=qi;
                vzB(j)=0;
                vregime(j)=1;
            elseif delta>=delta23
                vq2(j)=qstar;
                vy(j)=yi;
                vz2(j)=qstar/(1+i);
                vzB(j)=yi-(1-delta)*qstar/(1+i);
                vregime(j)=3;
            else
                x0=ystar;
                [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
                vq2(j)=q2;
                vy(j)=y;
                vz2(j)=y/(1-delta);
                vzB(j)=0;
                x0=ystar;
                vregime(j)=2;
            end
        end
        
        
        vzA=(1-vdelta).*vz2;
        vz=vzA+vzB;
        vzT=(1-vdelta).*vz2+vy;
        vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+vy+3*Xstar; % count CM
        
        vvelocity=vzT./vz;
        %     vtheta=vzT./vGDP;
        omega=0.0;
        vtheta=vzT./(vGDP-(1-omega)*3*Xstar);
        vrho=vz./vGDP;
        
    elseif model==3
        if flag==2
            lb=[0.0001,0.1,0,0.0001];
            x0=[0.5,2,10.0,0.5];
        elseif flag==3
            lb=[0.0001,0.2,0,0.0001];
            x0=[0.5,1.5,20,3]; %[eta s Xstar alpha]
        elseif flag==0
            x0=[0.5,2,10.0,0.8];
            lb=[0.0001,0.01,0,0.0001];
        else
            x0=[0.5,2,10.0,0.8];
            lb=[0.0001,0.1,0,0.0001];
        end
        A=[]; b=[]; Aeq=[]; beq=[];
        %     lb=[0.001,1,0,0.999];
        %     ub=[0.9999999,inf,inf,0.9999999];
        ub=[10,inf,inf,10.0];
        [x,ess,exitflag] = fmincon(@eqn_noconnfit1,x0,A,b,Aeq,beq,lb,ub,nonlcon,fmcopt);
        eta=x(1); s=x(2); Xstar=x(3);  alpha=x(4);
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
        
        for j=1:length(year)
            i=ir(j)/100;
            qi=(s/(1+i))^(1/alpha);
            yi=(1+i)^(-1/eta);
            yi3(j)=yi;
            qi3(j)=qi;
            qstar=s^(1/alpha);
            delta=vdelta(j);
            if delta>1
                vdelta(j)=1;
                delta=1;
            end
            vz(j)=(1-delta)*qi+yi;
            vGDP(j)=vz(j)+delta*qstar+4*Xstar;
        end
        vzT=vz;
        omega=0.0;
        vtheta=vzT./(vGDP-(1-omega)*4*Xstar);
        vrho=vz./vGDP;
        vtheta1=(vzT-yi3')./(vGDP-(1-omega)*4*Xstar-yi3');
    elseif model==4
        if flag==3
            x0=[0.5,1.5,100]; %[eta s Xstar]
        else
            x0=[0.8,2,100]; %[eta s Xstar]
        end
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0.001,0.001,0];
        ub=[1,inf,inf];
        [x,ess,exitflag] = fmincon(@eqn_noconnfit,x0,A,b,Aeq,beq,lb,ub,nonlcon,fmcopt);
        eta=x(1); s=x(2); Xstar=x(3);  %alpha=x(4);
        alpha=eta;
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
        
        for j=1:length(year)
            i=ir(j)/100;
            qi=(s/(1+i))^(1/alpha);
            yi=(1+i)^(-1/eta);
            yi3(j)=yi;
            qi3(j)=qi;
            qstar=s^(1/alpha);
            delta=vdelta(j);
            if delta>1
                vdelta(j)=1;
                delta=1;
            end
            vz(j)=(1-delta)*qi+yi;
            vGDP(j)=vz(j)+delta*qstar+4*Xstar;
        end
        vzT=vz;
        omega=0.0;
        vtheta=vzT./(vGDP-(1-omega)*4*Xstar);
        vrho=vz./vGDP;
        vtheta1=(vzT-yi3')./(vGDP-(1-omega)*4*Xstar-yi3');
        
    elseif model==5
        if flag==3
            x0=[0.9 2 10]; % [eta s Xstar]
        else
            x0=[0.9,1.5,5];
        end
        % x0=[0.9,1.5,5,0.5];
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0,1.01,0];
        ub=[0.999999,inf,inf];
        % lb=[0.01,1,0,0.99];
        % ub=[0.999999,inf,inf,0.999999];
        
        [x,ess,exitflag] = fmincon(@eqn_tfit2,x0,A,b,Aeq,beq,lb,ub,nonlcon,fmcopt);
        eta=x(1); s=x(2); Xstar=x(3); alpha=eta;
        %     alpha=x(4);
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
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
            delta12=1-ystar/qi;
            delta23=1-yi*(1+i)/qstar;
            
            if delta12<0
                NoRegime1=1;
            end
            if delta<=delta12
                vq2(j)=qi;
                vy(j)=ystar;
                vz2(j)=qi;
                vzB(j)=0;
                vregime(j)=1;
            elseif delta>=delta23
                vq2(j)=qstar;
                vy(j)=yi;
                vz2(j)=qstar/(1+i);
                vzB(j)=yi-(1-delta)*qstar/(1+i);
                vregime(j)=3;
            else
                x0=ystar;
                [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
                vq2(j)=q2;
                vy(j)=y;
                vz2(j)=y/(1-delta);
                vzB(j)=0;
                x0=ystar;
                vregime(j)=2;
            end
        end
        
        vzA=(1-vdelta).*vz2;
        vz=vzA+vzB;
        vzT=(1-vdelta).*vz2;
        vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+3*Xstar; % count CM
        vrho=vz./vGDP;
        vtheta=vzT./(vGDP-3*Xstar);
    elseif model==6
        if flag==3
            x0=[0.5,1.5,100]; %[eta s Xstar]
        else
            x0=[0.8,2,100]; %[eta s Xstar]
        end
        A=[]; b=[]; Aeq=[]; beq=[];
        lb=[0.001,0.001,0];
        ub=[1,inf,inf];
        [x,ess,exitflag] = fmincon(@eqn_noconnfit2,x0,A,b,Aeq,beq,lb,ub,nonlcon,fmcopt);
        eta=x(1); s=x(2); Xstar=x(3);  %alpha=x(4);
        alpha=eta;
        rmsd=sqrt(ess/n);    % root of mean squared error
        rmsd=rmsd/mean(rho);  % normalized RMSE
        
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
        
    end
    
    % Generate graphs
    
    f1=figure(1);
    subplot(1,2,1)
    plot(yt,ddata,'o');
    hold on
    plot(dt,'predobs');
    ylabel('\delta')
    xlabel('year')
    hold off
    
    subplot(1,2,2)
    [ir_t,gof1]=fit(year,ir,'poly2','Normalize','on','Robust','LAR');
    plot(year,ir,'o');
    hold on
    plot(ir_t,'predobs');
    ylabel('N. Int. Rate')
    xlabel('year')
    hold off
    [rhot,gof2]=fit(year,rho,'poly2','Normalize','on','Robust','LAR');
    [vrhot,gof3]=fit(year,vrho,'poly2','Normalize','on','Robust','LAR');
    
    f2=figure(2);
    subplot(2,2,1)
    p1=plot(year,rho,'o',year,vrho,'*r');
    hold on
    p2=plot(rhot,'b');
    p3=plot(vrhot,'-.r');
    ylabel('\rho')
    xlabel('year')
    hold off
    legend([p2,p3],'Data Trend','Model Trend','Location','best');
    ylim([0.025,0.075]);
    if (flag==2)
        xlim([1965 2015]);
    elseif (flag==3)
        xlim([1970 2015]);
        
    else
        xlim([1960 2015]);
    end
    
    [rhot,gof4]=fit(rho,vrho,'poly1')
    cor_rho=corr(rho,vrho);
    subplot(2,2,2)
    p1=plot(rho,vrho,'.b');
    hold on
    p2=plot(rhot);
    ylabel('\rho: Model')
    xlabel('\rho: Data')
    hold off
    legend(p2,'Regression Line','Location','best');
    
    [vthetat,gof5]=fit(year,vtheta,'poly2','Normalize','on','Robust','LAR');
    
    subplot(2,2,3)
    p1=plot(year,vtheta,'*r');
    hold on
    p2=plot(vthetat,'r','predobs');
    ylabel('\theta');
    xlabel('year');
    ylim([0 1])
    if (flag==2)
        xlim([1965 2015]);
    elseif (flag==3)
        xlim([1970 2015]);
        
    else
        xlim([1960 2015]);
    end
    
    [trt1,gof6]=fit(vrho(1:sn-1),vtheta(1:sn-1),'poly1');
    [trt2,gof6]=fit(vrho(sn:n),vtheta(sn:n),'poly1');
    trcorr1=corrcoef(vrho(1:sn-1),vtheta(1:sn-1));
    trcorr2=corrcoef(vrho(sn:n),vtheta(sn:n));
    lg1=strcat('Corr.Coef. before 1990=',num2str(trcorr1(1,2),2));
    lg2=strcat('Corr.Coef. after 1990=',num2str(trcorr2(1,2),2));
    subplot(2,2,4)
    p1=plot(vrho(1:sn-1),vtheta(1:sn-1),'ob');
    hold on
    p2=plot(vrho(sn:n),vtheta(sn:n),'*r');
    xlabel('Simulated \rho')
    ylabel('Simulated \theta')
    hold off
    legend([p1,p2],lg1,lg2,'Location','best');
    
elseif fig==1
    
    for mp=0:2
        if mp==0
            model=0;
        elseif mp==1
            model=4;
        else
            model=1;
        end
        
        if irflag==1
            load data;
            if flag==0 % Canada
                year=can(:,1);
                ir=can(:,2);
                ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
                yt=[1960,1980,1999,2005,2012]';
                rho=can(:,3);
            elseif flag==1 % US
                year=usa(:,1);
                ir=usa(:,2);
                ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
                yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
                rho=usa(:,3);
            elseif flag==2 % Australia
                year=aus(:,1);
                ir=aus(:,2);
                ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                    0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                    0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                    0.90064194]';
                yt=[1965,1994:2011]';
                rho=aus(:,3);
            else  % UK
                uk(1:4,:)=[];
                year=uk(:,1);
                ir=uk(:,2);
                ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
                yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
                rho=uk(:,3);
            end
        else
            load ir
            if flag==0 % Canada
                yr0=ir_can(:,1);
                if irflag==2
                    ir0=ir_can(:,3);    % DLM of actual inflation
                else
                    ir0=ir_can(:,4);     % Term structure based
                end
                rho0=ir_can(:,6);
                idx=(isnan(ir0))|(isnan(rho0));
                year=yr0(~idx);
                ir=ir0(~idx)*100;
                rho=rho0(~idx);
                ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
                yt=[1960,1980,1999,2005,2012]';
            elseif flag==1 % US
                yr0=ir_usa(:,1);
                if irflag==2
                    ir0=ir_usa(:,3);    % DLM of actual inflation
                else
                    ir0=ir_usa(:,4);    % Term structure based
                end
                rho0=ir_usa(:,7);
                idx=(isnan(ir0))|(isnan(rho0));
                year=yr0(~idx);
                ir=ir0(~idx)*100;
                rho=rho0(~idx);
                ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
                yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
            elseif flag==2 % Australia
                yr0=ir_aus(:,1);
                if irflag==2
                    ir0=ir_aus(:,3);    % DLM of actual inflation
                else
                    ir0=ir_aus(:,4);    % Term structure based
                end
                rho0=ir_aus(:,7);
                idx=(isnan(ir0))|(isnan(rho0));
                year=yr0(~idx);
                ir=ir0(~idx)*100;
                rho=rho0(~idx);
                ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                    0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                    0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                    0.90064194]';
                yt=[1960,1994:2011]';
            else  % UK
                yr0=ir_uk(:,1);
                if irflag==2
                    ir0=ir_uk(:,3);    % DLM of actual inflation
                else
                    ir0=ir_uk(:,4);    % Term structure based
                end
                rho0=ir_uk(:,6);
                idx=(isnan(ir0))|(isnan(rho0));
                year=yr0(~idx);
                ir=ir0(~idx)*100;
                rho=rho0(~idx);
                year(1:4)=[];
                ir(1:4)=[];
                rho(1:4)=[];
                ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
                yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
            end
        end
        % Use the data to find the fitted curve of delta
        [dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
        vdelta=feval(dt,year); % evaluate delta at each year between 1960 to 2014
        n=length(year);
        vq2=zeros(n,1);
        vy=vq2;
        vz2=vq2;
        vz=vz2;
        vzB=vq2;
        vregime=vq2;
        vGDP=vq2;
        vrho=vq2;
        vtheta=vq2;
        vzT=vq2;
        
        options=optimset('MaxFunEvals',50000 ,'Display', 'off', 'MaxIter',5000);
        nonlcon=[];
        fmcopt=optimoptions(@fmincon,'MaxFunEvals',50000,'MaxIter',5000);% Option to display output
        if model==1
            if flag==3
                x0=[0.9 2 10]; % [eta s Xstar] initial guess
            else
                x0=[0.9,1.5,5];
            end
            A=[]; b=[]; Aeq=[]; beq=[];
            lb=[0,1.01,0];
            ub=[0.999999,inf,inf];
            
            [x,ess,exitflag] = fmincon(@eqn_tfit,x0,A,b,Aeq,beq,lb,ub);
            eta=x(1); s=x(2); Xstar=x(3); alpha=eta;
            rmsd=sqrt(ess/n);    % root of mean squared error
            rmsd=rmsd/mean(rho);  % normalized RMSE
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
                delta12=1-ystar/qi; % Cutoff of Regimes 1 and 2
                delta23=1-yi*(1+i)/qstar; % Cutoff of Regimes 2 and 3
                
                if delta12<0
                    NoRegime1=1;
                end
                if delta<=delta12   % Regime 1
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
                else                % Regime 2
                    x0=ystar;
                    [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
                    vq2(j)=q2;
                    vy(j)=y;
                    vz2(j)=y/(1-delta);
                    vzB(j)=0;
                    x0=ystar;
                    vregime(j)=2;
                end
            end
            
            
            vzA=(1-vdelta).*vz2;
            vz=vzA+vzB;
            vzT=(1-vdelta).*vz2+vy;
            vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+vy+3*Xstar; % count CM
            
            vvelocity=vzT./vz;
            vtheta=vzT./(vGDP-3*Xstar);
            vrho=vz./vGDP;
            
        elseif model==0
            if flag==2
                x0=[0.5,1.5,15];
            else
                x0=[0.9,1.5,10.0];
            end
            A=[]; b=[]; Aeq=[]; beq=[];
            lb=[0,1,0];
            ub=[0.99999,1,inf];
            [x,ess,exitflag] = fmincon(@eqn_LWfit,x0,A,b,Aeq,beq,lb,ub);
            rmsd=sqrt(ess/n);    % root of mean squared error
            rmsd=rmsd/mean(rho);  % normalized RMSE
            
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
            
        elseif model==4
            if flag==3
                x0=[0.5,1.5,100]; %[eta s Xstar]
            else
                x0=[0.8,2,100]; %[eta s Xstar]
            end
            A=[]; b=[]; Aeq=[]; beq=[];
            lb=[0.001,0.001,0];
            ub=[1,inf,inf];
            [x,ess,exitflag] = fmincon(@eqn_noconnfit,x0,A,b,Aeq,beq,lb,ub,nonlcon,fmcopt);
            eta=x(1); s=x(2); Xstar=x(3);  %alpha=x(4);
            alpha=eta;
            rmsd=sqrt(ess/n);    % root of mean squared error
            rmsd=rmsd/mean(rho);  % normalized RMSE
            
            for j=1:length(year)
                i=ir(j)/100;
                qi=(s/(1+i))^(1/alpha);
                yi=(1+i)^(-1/eta);
                yi3(j)=yi;
                qi3(j)=qi;
                qstar=s^(1/alpha);
                delta=vdelta(j);
                if delta>1
                    vdelta(j)=1;
                    delta=1;
                end
                vz(j)=(1-delta)*qi+yi;
                vGDP(j)=vz(j)+delta*qstar+4*Xstar;
            end
            vzT=vz;
            omega=0.0;
            vtheta=vzT./(vGDP-(1-omega)*4*Xstar);
            vrho=vz./vGDP;
            vtheta1=(vzT-yi3')./(vGDP-(1-omega)*4*Xstar-yi3');
            
        end
        
        f1=figure(1);
        subplot(2,3,mp+1)
        [rhot,gof2]=fit(year,rho,'poly2','Normalize','on','Robust','LAR');
        [vrhot,gof3]=fit(year,vrho,'poly2','Normalize','on','Robust','LAR');
        p1=plot(year,rho,'o',year,vrho,'*r');
        hold on
        p2=plot(rhot,'b');
        p3=plot(vrhot,'--r');
        ylabel('\rho')
        xlabel('year')
        hold off
        legend([p2,p3],'Data Trend','Model Trend','Location','best');
        %ylim([0.025,0.075]);
        if (flag==2)
            xlim([1965 2015]);
        elseif (flag==3)
            xlim([1970 2015]);
            
        else
            xlim([1960 2015]);
        end
        
        [vthetat,gof5]=fit(year,vtheta,'poly2','Normalize','on','Robust','LAR');
        
        subplot(2,3,mp+4)
        p1=plot(year,vtheta,'*r');
        hold on
        p2=plot(vthetat,'r','predobs');
        ylabel('\theta');
        xlabel('year');
        ylim([0 1])
        if (flag==2)
            xlim([1965 2015]);
        elseif (flag==3)
            xlim([1970 2015]);
            
        else
            xlim([1960 2015]);
        end
        
        [trt1,gof6]=fit(vrho(1:sn-1),vtheta(1:sn-1),'poly1');
        [trt2,gof6]=fit(vrho(sn:n),vtheta(sn:n),'poly1');
        trcorr1=corr(vrho(1:sn-1),vtheta(1:sn-1));
        trcorr2=corr(vrho(sn:n),vtheta(sn:n));
        lg1=strcat('Corr.Coef. before 1990=',num2str(trcorr1,2));
        lg2=strcat('Corr.Coef. after 1990=',num2str(trcorr2,2));
        
        f2=figure(2);
        [rhot,gof4]=fit(rho,vrho,'poly1')
        cor_rho=corr(rho,vrho);
        subplot(2,3,mp+1)
        p1=plot(rho,vrho,'.b');
        hold on
        p2=plot(rhot);
        ylabel('\rho: Model')
        xlabel('\rho: Data')
        hold off
        legend(p2,'Regression Line','Location','best');
        
        subplot(2,3,mp+4)
        p1=plot(vrho(1:sn-1),vtheta(1:sn-1),'ob');
        hold on
        p2=plot(vrho(sn:n),vtheta(sn:n),'*r');
        xlabel('Simulated \rho')
        ylabel('Simulated \theta')
        ylim([0 1])
        hold off
        legend([p1,p2],lg1,lg2,'Location','best');
    end
else
    for irflag=1:3
        for flag=0:3
            if irflag==1
                load data;
                if flag==0 % Canada
                    year=can(:,1);
                    ir=can(:,2);
                    ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
                    yt=[1960,1980,1999,2005,2012]';
                    rho=can(:,3);
                elseif flag==1 % US
                    year=usa(:,1);
                    ir=usa(:,2);
                    ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
                    yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
                    rho=usa(:,3);
                elseif flag==2 % Australia
                    year=aus(:,1);
                    ir=aus(:,2);
                    ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                        0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                        0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                        0.90064194]';
                    yt=[1965,1994:2011]';
                    rho=aus(:,3);
                else  % UK
                    uk(1:4,:)=[];
                    year=uk(:,1);
                    ir=uk(:,2);
                    ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
                    yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
                    rho=uk(:,3);
                end
            else
                load ir
                if flag==0 % Canada
                    yr0=ir_can(:,1);
                    if irflag==2
                        ir0=ir_can(:,3);    % DLM of actual inflation
                    else
                        ir0=ir_can(:,4);     % Term structure based
                    end
                    rho0=ir_can(:,6);
                    idx=(isnan(ir0))|(isnan(rho0));
                    year=yr0(~idx);
                    ir=ir0(~idx)*100;
                    rho=rho0(~idx);
                    ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
                    yt=[1960,1980,1999,2005,2012]';
                elseif flag==1 % US
                    yr0=ir_usa(:,1);
                    if irflag==2
                        ir0=ir_usa(:,3);    % DLM of actual inflation
                    else
                        ir0=ir_usa(:,4);    % Term structure based
                    end
                    rho0=ir_usa(:,7);
                    idx=(isnan(ir0))|(isnan(rho0));
                    year=yr0(~idx);
                    ir=ir0(~idx)*100;
                    rho=rho0(~idx);
                    ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
                    yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
                elseif flag==2 % Australia
                    yr0=ir_aus(:,1);
                    if irflag==2
                        ir0=ir_aus(:,3);    % DLM of actual inflation
                    else
                        ir0=ir_aus(:,4);    % Term structure based
                    end
                    rho0=ir_aus(:,7);
                    idx=(isnan(ir0))|(isnan(rho0));
                    year=yr0(~idx);
                    ir=ir0(~idx)*100;
                    rho=rho0(~idx);
                    ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
                        0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
                        0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
                        0.90064194]';
                    yt=[1960,1994:2011]';
                else  % UK
                    yr0=ir_uk(:,1);
                    if irflag==2
                        ir0=ir_uk(:,3);    % DLM of actual inflation
                    else
                        ir0=ir_uk(:,4);    % Term structure based
                    end
                    rho0=ir_uk(:,6);
                    idx=(isnan(ir0))|(isnan(rho0));
                    year=yr0(~idx);
                    ir=ir0(~idx)*100;
                    rho=rho0(~idx);
                    year(1:4)=[];
                    ir(1:4)=[];
                    rho(1:4)=[];
                    ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
                    yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
                end
            end
            % Use the data to find the fitted curve of delta
            [dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
            vdelta=feval(dt,year); % evaluate delta at each year between 1960 to 2014
            n=length(year);
            vq2=zeros(n,1);
            vy=vq2;
            vz2=vq2;
            vz=vz2;
            vzB=vq2;
            vregime=vq2;
            vGDP=vq2;
            vrho=vq2;
            vtheta=vq2;
            vzT=vq2;
            
            options=optimset('MaxFunEvals',50000 ,'Display', 'off', 'MaxIter',5000);
            nonlcon=[];
            fmcopt=optimoptions(@fmincon,'MaxFunEvals',50000,'MaxIter',5000);% Option to display output
            if flag==3
                x0=[0.9 2 10]; % [eta s Xstar] initial guess
            else
                x0=[0.9,1.5,5];
            end
            A=[]; b=[]; Aeq=[]; beq=[];
            lb=[0,1.01,0];
            ub=[0.999999,inf,inf];
            
            [x,ess,exitflag] = fmincon(@eqn_tfit,x0,A,b,Aeq,beq,lb,ub);
            eta=x(1); s=x(2); Xstar=x(3); alpha=eta;
            rmsd=sqrt(ess/n);    % root of mean squared error
            rmsd=rmsd/mean(rho);  % normalized RMSE
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
                delta12=1-ystar/qi; % Cutoff of Regimes 1 and 2
                delta23=1-yi*(1+i)/qstar; % Cutoff of Regimes 2 and 3
                
                if delta12<0
                    NoRegime1=1;
                end
                if delta<=delta12   % Regime 1
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
                else                % Regime 2
                    x0=ystar;
                    [x,fval,exitflag] = fzero(@eqn_Regime201610,x0,options);
                    vq2(j)=q2;
                    vy(j)=y;
                    vz2(j)=y/(1-delta);
                    vzB(j)=0;
                    x0=ystar;
                    vregime(j)=2;
                end
            end
            
            
            vzA=(1-vdelta).*vz2;
            vz=vzA+vzB;
            vzT=(1-vdelta).*vz2+vy;
            vGDP=((1-vdelta).*vq2+vdelta*qstar)./(vy.^(-eta))+vy+3*Xstar; % count CM
            
            vvelocity=vzT./vz;
            vtheta=vzT./(vGDP-3*Xstar);
            vrho=vz./vGDP;
            
            f1=figure(1);
            if flag==2
                mp=0;
            elseif flag==0
                mp=1;
            elseif flag==3
                mp=2;
            else
                mp=3;
            end
            subplot(4,3,3*mp+irflag)
            [rhot,gof2]=fit(year,rho,'poly2','Normalize','on','Robust','LAR');
            [vrhot,gof3]=fit(year,vrho,'poly2','Normalize','on','Robust','LAR');
            p1=plot(year,rho,'o',year,vrho,'*r');
            hold on
            p2=plot(rhot,'b');
            p3=plot(vrhot,'--r');
            ylabel('\rho')
            xlabel('year')
            hold off
            legend([p2,p3],'Data Trend','Model Trend','Location','best');
            %ylim([0.025,0.075]);
            if (flag==2)
                xlim([1965 2015]);
            elseif (flag==3)
                xlim([1970 2015]);
                
            else
                xlim([1960 2015]);
            end
            if irflag==1
                
                f4=figure(4);
                if flag==2
                    ax=subplot(2,2,1);
                    plot(ax,year,vregime,'*r');
                    title('Australia');
                elseif flag==0
                    ax=subplot(2,2,2);
                    plot(ax,year,vregime,'*r');
                    title('Canada');
                elseif flag==3
                    ax=subplot(2,2,3);
                    plot(ax,year,vregime,'*r');
                    title('UK');
                else
                    ax=subplot(2,2,4);
                    plot(ax,year,vregime,'*r');
                    title('US');
                end
                ylabel('Regime')
                xlabel('year')
                ylim([1 3]);
                ax.YTick=[1,2,3];
                
            end
            
        end
    end
end

