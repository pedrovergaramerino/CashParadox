module CashParadox
export load_data
using MAT
using Statistics
using LinRegOutliers
using Plots
using Econometrics
using Roots

"""
   load_data(irflag,flag)

Load data depending on the figure and country wanted.

"""
mutable struct load_data
    year::Vector{Float64}
    ir::Vector{Float64}
    yt::Vector{Int64}
    ρ::Vector{Float64}
    ddata::Vector{Float64}

function load_data(irflag::Int64,flag::Int64)
    dict=matread("JiangShaoCodeData\\data.mat")
    dictir=matread("JiangShaoCodeData\\ir.mat")
    can=dict["can"]
    uk=dict["uk"]
    usa=dict["usa"]
    aus=dict["aus"]
    ir_can=dictir["ir_can"]
    ir_uk=dictir["ir_uk"]
    ir_usa=dictir["ir_usa"]
    ir_aus=dictir["ir_aus"]
    if irflag==1
        if flag==0 #CANADA
            year=can[:,1]
            ir=can[:,2]
            ddata=[0 0.49 0.782 0.822 0.856]' #Data on access to credit
            yt=[1960 1980 1999 2005 2012]'
            ρ=can[:,3]
        elseif flag==1 #USA
            year=usa[:,1]
            ir=usa[:,2]
            ddata=[0 16.3 38.3 43.0 55.8 62.2 66.5 68 72.6 71.5]'/100
            yt=[1960 1970 1977 1983 1989 1992 1995 1998 2001 2004]'
            ρ=usa[:,3]
        elseif flag==2 #AUSTRALIA
            year=aus[:,1]
            ir=aus[:,2]
            ddata=[0 0.520354595 0.537666422 0.57032816 0.594028535 0.622316221 0.65072652 0.687205136 0.694317256 0.73605656 0.760156768 0.792611398 0.841961623 0.878037024 0.905946326 0.90592637 0.897747613 0.904361411 0.90064194]'
            yt= vcat([1965], collect(1994:2011))
            ρ=aus[:,3]
        else #UK 
            uk=uk[5:end,:]
            year=uk[:,1]
            ir=uk[:,2]
            ddata=[0 0.5 0.53 0.65 0.62 0.64 0.62 0.61  0.6]'
            yt=[1965 1999  2000  2005  2009  2010  2011  2012  2013]'
            ρ=uk[:,3]
        end
    else
        if flag==0 #CANADA
            yr0=ir_can[:,1]
            if irflag==2
                ir0=ir_can[:,3] #DLM of actual inflation
            else
                ir0=ir_can[:,4] #Term structure based
            end
            ρ0=ir_can[:,6]
            idx=Vector{Bool}(undef,length(ir0))
            for i in 1:length(ir0)
                if (isnan(ρ0[i])|isnan(ir0[i]))
                    idx[i]=0
                else
                    idx[i]=1
                end
            end
            year=yr0[idx]
            ir=ir0[idx]*100
            ρ=ρ0[idx]
            ddata=ddata=[0,0.49,0.782,0.822,0.856]' #Data on access to credit
            yt=[1960,1980,1999,2005,2012]'
        elseif flag==1 #USA
            yr0=ir_usa[:,1]
            if irflag==2
                ir0=ir_usa[:,3] #DLM of actual inflation
            else
                ir0=ir_usa[:,4] #Term structure based
            end
            ρ0=ir_usa[:,7]
            idx=Vector{Bool}(undef,length(ir0))
            for i in 1:length(ir0)
                if (isnan(ρ0[i])|isnan(ir0[i]))
                    idx[i]=0
                else
                    idx[i]=1
                end
            end
            year=yr0[idx]
            ir=ir0[idx]*100
            ρ=ρ0[idx]
            ddata=[0 16.3 38.3 43.0 55.8 62.2 66.5 68 72.6 71.5]'/100
            yt=[1960 1970 1977 1983 1989 1992 1995 1998 2001 2004]'
        elseif flag==2 #AUSTRALIA
            yr0=ir_aus[:,1]
            if irflag==2
                ir0=ir_aus[:,3] #DLM of actual inflation
            else
                ir0=ir_aus[:,4] #Term structure based
            end
            ρ0=ir_aus[:,7]
            idx=Vector{Bool}(undef,length(ir0))
            for i in 1:length(ir0)
                if (isnan(ρ0[i])|isnan(ir0[i]))
                    idx[i]=0
                else
                    idx[i]=1
                end
            end
            year=yr0[idx]
            ir=ir0[idx]*100
            ρ=ρ0[idx] 
            ddata=[0 0.520354595 0.537666422 0.57032816 0.594028535 0.622316221 0.65072652 0.687205136 0.694317256 0.73605656 0.760156768 0.792611398 0.841961623 0.878037024 0.905946326 0.90592637 0.897747613 0.904361411 0.90064194]'
            yt=vcat([1965], collect(1994:2011))  
        else #UK
            yr0=ir_uk[:,1]
            if irflag==2
                ir0=ir_uk[:,3]
            else
                ir0=ir_uk[:,4]
            end
            ρ0=ir_uk[:,6]
            for i in 1:length(ir0)
                if (isnan(ρ0[i])|isnan(ir0[i]))
                    idx[i]=0
                else
                    idx[i]=1
                end
            end
            year=yr0[idx]
            ir=ir0[idx]*100
            ρ=ρ0[idx]
            year=year[5:end]
            ir=ir[5:end]
            ρ=ρ[5:end]
            ddata=[0 0.5 0.53 0.65 0.62 0.64 0.62 0.61  0.6]'
            yt=[1965 1999  2000  2005  2009  2010  2011  2012  2013]'
        end
    end
    ddata=Vector{Float64}(vec(ddata))
    yt=Vector{Int64}(vec(yt))
    this=new()
    this.year=year
    this.ir=ir
    this.ρ=ρ
    this.ddata=ddata
    this.yt=yt
    return this
end
end


"""
   vδ(irflag,flag)

Computes vδ

"""
function vδ(irflag::Int64,flag::Int64)
    data=load_data(irflag,flag)
    normyt=(data.yt.- mean(data.yt))./std(data.yt)
    normyear=(data.year.- mean(data.yt))./std(data.yt)
    Xyt=[ones(length(normyt)) normyt normyt.^2]
    Xyear=[ones(length(normyear)) normyear normyear.^2]
    dt=lad(Xyt,data.ddata;starting_betas=nothing)
    vδ=Xyear*dt["betas"]
return vδ
end

"""
eqn_Regime201610

Extensive margin
Equations to solve q2 and y in regime 2

"""
function eqn_Regime201610(x::Float64,α::Float64,δ::Float64,η::Float64,s::Float64, i::Float64)
y=x
u_y=y^(-η)
q2=y*u_y/(1-δ)
v_q=s*q2^(-α)
F=u_y*v_q-(1+i)
return F
end

"""
eqn_tfit

Matching the time series of CIC_GDP ratio in the data with the simulated 
series in the full model

"""
function eqn_tfit(x::Vector{Float64}, irflag::Int64, flag::Int64)
    η=x[1];
    s=x[2];
    Xstar=x[3];
    α=η;
    data=load_data(irflag,flag)
    qstar=s^(1/α)
    n=length(data.year);
    ystar=1.0;
    vq2=zeros(n,1);
    vy=zeros(n,1);
    vz2=zeros(n,1);
    v_q=zeros(n,1)
    vz=zeros(n,1);
    vzT=zeros(n,1);
    vρ=zeros(n,1);
    vθ=zeros(n,1);
    vzA=zeros(n,1)
    vzB=zeros(n,1);
    vregime=zeros(n,1);
    vGDP=zeros(n,1);
    vδprime=vδ(irflag, flag)
    i=zeros(n,1)
    qi=zeros(n,1)
    yi=zeros(n,1)
    δ12=zeros(n,1)
    δ23=zeros(n,1)
    NoRegime1=zeros(n,1)
    y=zeros(n,1)
    u_y=zeros(n,1)
    q2=zeros(n,1)
    for j in 1:length(data.year)
        i[j]=data.ir[j]./100;
        qi[j]=(s./(1.0 .+i[j])).^(1.0 ./α);
        yi[j]=(1.0 .+i[j]) .^(-1.0 ./η);
        vδprime[j];
        if vδprime[j]>1
            vδprime[j]=1;
        end
        δ12[j]=1.0 .-ystar./qi[j];     # Cutoff of Regimes 1 and 2
        δ23[j]=1.0 .-yi[j] .*(1.0 .+i[j])./qstar;   #Cutoff of Regimes 2 and 3
        
        if δ12[j]<0.0
            NoRegime1[j]=1;
        end
        if vδprime[j]<=δ12[j]       # Regime 1
            vq2[j]=qi[j];
            vy[j]=ystar;
            vz2[j]=qi[j];
            vzB[j]=0.0;
            vregime[j]=1;
        elseif vδprime[j]>=δ23[j]   # Regime 3
            vq2[j]=qstar;
            vy[j]=yi[j];
            vz2[j]=qstar./(1.0 .+i[j]);
            vzB[j]=yi[j] .-(1.0 .-vδprime[j]) .*qstar./(1.0 .+i[j]);
            vregime[j]=3;
        else                    # Regime 2
            y[j] = find_zero( x-> eqn_Regime201610(x,α,vδprime[j],η,s,i[j]),ystar);
            u_y[j]=y[j].^(-η);
            q2[j]=y[j].*u_y[j]./(1.0 .-vδprime[j]);
            v_q[j]=s.*q2[j].^(-α)
            vq2[j]=q2[j];
            vy[j]=y[j];
            vz2[j]=y[j]./(1.0 .-vδprime[j]);
            vzB[j]=0.0;
            vregime[j]=2;
        end

    end
    vzA.=(1 .-vδprime).*vz2;
    vz.=vzA .+vzB;
    vzT.=(1 .-vδprime).*vz2+vy;
    vGDP.=((1 .-vδprime).*vq2+vδprime.*qstar)./(vy.^(-η)).+vy.+3 .*Xstar; # count CM
    vρ.=vz./vGDP;
    ess=sum((vρ .-data.ρ).^2);
    return ess

end

"""
eqn_LWfit

Matching the time series of CIC_GDP ratio in the data with the simulated 
series in the full model

"""
function eqn_LWfit(x::Vector{Float64}, irflag::Int64, flag::Int64)
    η=x[1];
    s=x[2];
    Xstar=x[3];
    α=η;
    data=load_data(irflag,flag)
    qstar=s^(1/α)
    n=length(data.year);
    ystar=1.0;
    vz=zeros(n,1);
    vzT=zeros(n,1);
    vρ=zeros(n,1);
    vθ=zeros(n,1);
    vGDP=zeros(n,1);
    vδprime=vδ(irflag, flag)
    i=zeros(n,1)
    qi=zeros(n,1)
    for j in 1:length(data.year)
        i[j]=data.ir[j]./100;
        qi[j]=(s./(1.0 .+i[j])).^(1.0 ./α);
        if  vδprime[j]>1
            vδprime[j]=1;
        end
        vz[j]=(1.0 .- vδprime[j]).* qi[j]
        vGDP[j]=(1.0 .- vδprime[j]).* qi[j] + vδprime[j].*qstar+2.0 .*Xstar
    end

    vzT=vz
    ω=0.0
    vθ=vzT./(vGDP.-(1.0 .- ω).*2.0.*Xstar)
    vρ=vz./vGDP
    ess=sum((vρ.-data.ρ).^2)
    return ess

end

"""
eqn_noconnfit

Matching the time series of CIC_GDP ratio in the data with the simulated 
series in the full model

"""
function eqn_noconnfit(x::Vector{Float64}, irflag::Int64, flag::Int64)
    η=x[1];
    s=x[2];
    Xstar=x[3];
    α=η;
    data=load_data(irflag,flag)
    qstar=s^(1/α)
    n=length(data.year);
    vz=zeros(n,1);
    vzT=zeros(n,1);
    vρ=zeros(n,1);
    vθ=zeros(n,1);
    vGDP=zeros(n,1);
    vδprime=vδ(irflag, flag)
    i=zeros(n,1)
    qi=zeros(n,1)
    yi=zeros(n,1)
    for j in 1:length(data.year)
        i[j]=data.ir[j]./100;
        qi[j]=(s./(1.0 .+i[j])).^(1.0 ./α);
        yi[j]=(1.0.+i[j]).^(-1.0./η)
        if  vδprime[j]>1
            vδprime[j]=1;
        elseif vδprime[j]<0
            vδprime[j]=0
        end
        vz[j]=(1.0 .- vδprime[j]).* qi[j] .+yi[j]
        vGDP[j]=vz[j]+vδprime[j].*qstar+4.0*Xstar
    end

    vzT=vz
    ω=0.0
    vθ=vzT./(vGDP.-(1.0 .- ω).*4.0.*Xstar)
    vρ=vz./vGDP
    ess=sum((vρ.-data.ρ).^2)
    return ess

end

"""
vs(irflag,flag,model)

Returns vθ and vρ for a specified model

"""
function vs(irflag::Int64, flag::Int64, model::Int64)
    if model==1
        if flag==3
            x0=[0.9 2.0 10.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.9,1.5,5.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.0,1.01,0.0];
        ub=[0.999999,Inf,Inf];
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vq2=zeros(n,1);
        vy=zeros(n,1);
        vz2=zeros(n,1);
        vz=zeros(n,1);
        vzB=zeros(n,1);
        vregime=zeros(n,1);
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt= fmincon(x-> eqn_tfit(x,irflag,flag), x0,R,r,lb,ub; tol = 1e-15, iterlim=50000);
        x=opt[1]
        ess=opt[2]
        η=x[1]; 
        s=x[2]; 
        Xstar=x[3]; 
        α=η;
        qstar=s^(1/α)
        rmsd=sqrt(ess/n);    #root of mean squared error
        rmsd=rmsd/mean(data.ρ);  # normalized RMSE
        for j=1:length(data.year)
            i=data.ir[j]/100;
            qi=(s/(1+i))^(1/α);
            yi=(1+i)^(-1/η);
            δ=vδprime[j];
            if δ>1
                vδprime[j]=1;
                δ=1;
            end
            δ12=1-ystar/qi; # Cutoff of Regimes 1 and 2
            δ23=1-yi*(1+i)/qstar; # Cutoff of Regimes 2 and 3
            
            if δ12<0
                NoRegime1=1;
            end
            if δ<=δ12   # Regime 1
                vq2[j]=qi;
                vy[j]=ystar;
                vz2[j]=qi;
                vzB[j]=0;
                vregime[j]=1;
            elseif δ>=δ23   # Regime 3
                vq2[j]=qstar;
                vy[j]=yi;
                vz2[j]=qstar/(1+i);
                vzB[j]=yi-(1-δ)*qstar/(1+i);
                vregime[j]=3;
            else                # Regime 2
                x0=ystar;
                eqn_Regime2(x)=eqn_Regime201610(x,α,δ,η,s,i)
                y = find_zero(eqn_Regime2,x0);
                u_y=y^(-η);
                q2=y*u_y/(1-δ);
                v_q=s*q2^(-α)
                vq2[j]=q2;
                vy[j]=y;
                vz2[j]=y/(1-δ);
                vzB[j]=0;
                vregime[j]=2;
            end
        end
        
        vzA=(1.0 .-vδprime).*vz2;
        vz=vzA .+ vzB;
        vzT=(1.0 .-vδprime).*vz2 .+vy;
        vGDP=((1.0 .-vδprime).*vq2.+vδprime .*qstar)./(vy.^(-η)) .+vy .+3 .*Xstar; # count CM
        
        vvelocity=vzT./vz;
        vθ=vzT./(vGDP .-3 .*Xstar);
        vρ=vz./vGDP;
        return (vθ,vρ)
    
    elseif model==0
        if flag==2
            x0=[0.5 1.5 15.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.9,1.5,10.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.0,1.0,0.0];
        ub=[0.999999,1.0,Inf];
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt= fmincon(x-> eqn_tfit(x,irflag,flag), x0,R,r,lb,ub; tol = 1e-15, iterlim=50000);
        x=opt[1]
        ess=opt[2]
        η=x[1]; 
        s=x[2]; 
        Xstar=x[3]; 
        α=η;
        qstar=s^(1/α)
        rmsd=sqrt(ess/n);    #root of mean squared error
        rmsd=rmsd/mean(data.ρ);
        for j=1:length(data.year)
            i=data.ir[j]/100;
            qi=(s./(1+i))^(1.0./α);
            qstar=s^(1.0./α);
            if vδprime[j]>1
                vδprime[j]=1;
            end
            vz[j]=(1.0.-vδprime[j]).*qi;
            vGDP[j]=(1.0.-vδprime[j])*qi.+vδprime[j].*qstar.+2.0.*Xstar;
        end
        vzT=vz;
        ω=0.0;
        vθ=vzT./(vGDP.-(1.0.-ω).*2.0.*Xstar);
        vρ=vz./vGDP; 
        return (vθ,vρ)       

    elseif model==4
        if flag==3
            x0=[0.5 1.5 100.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.8,2.0,100.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.001,0.001,0.0];
        ub=[1.0,Inf,Inf];
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        vθ1=zeros(n,1)
        yi3=zeros(n,1)
        qi3=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt= fmincon(x-> eqn_noconnfit(x,irflag,flag), x0,R,r,lb,ub; tol = 1e-15, iterlim=50000);
        x=opt[1]
        ess=opt[2]
        η=x[1]; 
        s=x[2]; 
        Xstar=x[3]; 
        α=η;
        qstar=s^(1/α)
        rmsd=sqrt(ess/n);    #root of mean squared error
        rmsd=rmsd/mean(data.ρ);
        for j=1:length(data.year)
            i=data.ir[j]./100.0;
            qi=(s./(1.0.+i)).^(1.0./α);
            yi=(1.0.+i).^(-1.0./η);
            yi3[j]=yi;
            qi3[j]=qi;
            if vδprime[j]>1
               vδprime[j]=1;
            end
            vz[j]=(1.0.-vδprime[j]).*qi.+yi;
            vGDP[j]=vz[j].+vδprime[j].*qstar.+4.0.*Xstar;
        end
        vzT=vz;
        ω=0.0;
        vθ=vzT./(vGDP.-(1.0.-ω).*4.0.*Xstar);
        vρ=vz./vGDP;
        return (vθ,vρ)
    end  
end

"""
figsub5(irflag,flag)

Creates the subplots for figure 5 for a specific interest rate specification and country

"""
function figsub5(irflag::Int64 , flag::Int64, model::Int64)
    (vθ,vρ)=vs(irflag,flag,model)
    vθ=Vector{Float64}(vec(vθ))
    vρ=Vector{Float64}(vec(vρ))
    data=load_data(irflag,flag)
    normyear=(data.year.- mean(data.year))./std(data.year) 
    Xyear=[ones(length(normyear)) normyear normyear.^2]
    ladρt=lad(Xyear,data.ρ;starting_betas=nothing)
    betasρt=ladρt["betas"]
    fρt(x)=betasρt[1].+betasρt[2].*x.+betasρt[3].*x.^2
    ρt=[fρt(x) for x ∈ minimum(normyear):0.01:maximum(normyear)]
    ladvρt=lad(Xyear,vρ;starting_betas=nothing)
    betasvρt=ladvρt["betas"]
    fvρt(x)=betasvρt[1].+betasvρt[2].*x.+betasvρt[3].*x.^2
    vρt=[fvρt(x) for x ∈ minimum(normyear):0.01:maximum(normyear)]
    xgrid=minimum(normyear):0.01:maximum(normyear)
    xgrid=(xgrid.*std(data.year)).+mean(data.year)
    p=plot(xgrid,ρt,color=:lightskyblue,label="Data Trend",xlabel="year",ylabel="\\rho")
    scatter!(p,data.year,data.ρ,color=:lightskyblue,label=false)
    plot!(p,xgrid,vρt,color=:red,style=:dash,label="Model Trend")
    scatter!(p,data.year,vρ,color=:red,shape=:star5,label=false)
return p
end


end #module
