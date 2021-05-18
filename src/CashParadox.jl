module CashParadox
using MAT
using Statistics
using LinRegOutliers
using Plots
using Roots
using Optim
using DataFrames

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
    dict=matread("JiangShaoCodeData//data.mat")
    dictir=matread("JiangShaoCodeData//ir.mat")
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
            ddata=[0.0, 16.3, 38.3, 43.0, 55.8, 62.2, 66.5, 68.0, 72.6, 71.5]'/100
            yt=[1960, 1970, 1977, 1983, 1989, 1992, 1995, 1998 ,2001 ,2004]'
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
function eqn_Regime201610(x,α,δ,η,s,i)
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
function eqn_tfit(x, irflag::Int64, flag::Int64)
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
function eqn_noconnfit(x, irflag::Int64, flag::Int64)
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
    opt=[[0.0,0.0,0.0],0.0]
    if model==1
        if flag==3
            x0=[0.9 2.0 10.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        
        elseif (flag==1 && irflag==3)
            x0=[0.9,1.1,3.0]
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
        opt=Optim.optimize(x-> eqn_tfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        return (vθ,vρ,vregime)
    
    elseif model==0
        if flag==2
            x0=[0.5 1.0 15.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.9,1.0,10.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.0,0.9999999,0.0];
        ub=[0.999999,1.000001,Inf];
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vregime=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt=Optim.optimize(x-> eqn_LWfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        return (vθ,vρ,vregime)       

    elseif model==4
        if flag==3
            x0=[0.5 1.5 100.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.8,2.0,100.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.001,0.001,0.0];
        ub=[1.0,Inf,Inf]
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vregime=zeros(n,1)
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        yi3=zeros(n,1)
        qi3=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt=Optim.optimize(x-> eqn_noconnfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        return (vθ,vρ,vregime)
    end  
end

"""
figsub5(irflag,flag,model)

Creates the subplots for figure 5 for a specific interest rate specification and country

"""
function figsub5(irflag::Int64 , flag::Int64, model::Int64)
    (vθ,vρ,vregime)=vs(irflag,flag,model)
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
    p1=plot(xgrid,ρt,color=:lightskyblue,label="Data Trend",xlabel="year",ylabel="\\rho")
    scatter!(p1,data.year,data.ρ,color=:lightskyblue,label=false)
    plot!(p1,xgrid,vρt,color=:red,style=:dash,label="Model Trend")
    scatter!(p1,data.year,vρ,color=:red,shape=:cross,label=false)
    ladvθt=lad(Xyear,vθ;starting_betas=nothing)
    betasvθt=ladvθt["betas"]
    fvθt(x)=betasvθt[1].+betasvθt[2].*x.+betasvθt[3].*x.^2
    vθt=[fvθt(x) for x ∈ minimum(normyear):0.01:maximum(normyear)]
    p2=plot(xgrid,vθt,color=:red,label="fitted value",ylims=(0,1))
    scatter!(p2,data.year,vθ,color=:red,shape=:cross,label=false)
    p=plot(p1,p2,layout=(2,1))
return p
end

"""
fig5(flag)

Creates figure 5 for a specific country

"""
function fig5(flag::Int64)
p1=figsub5(1,flag,0)
p2=figsub5(1,flag,4)
p3=figsub5(1,flag,1)
p=plot(p1,p2,p3,layout=(1,3),size=(1000,1000))
return p
end

"""
figsubA2(irflag,flag,model)

Creates the subplots for figure 5 for a specific interest rate specification and country

"""
function figsubA2(irflag::Int64 , flag::Int64, model::Int64)
    if flag==2
        sn=20
    elseif flag==3
        sn=23
    else 
        sn=30
    end
    (vθ,vρ,vregime)=vs(irflag,flag,model)
    vθ=Vector{Float64}(vec(vθ))
    vρ=Vector{Float64}(vec(vρ))
    n=length(vθ)
    vθ1=vθ[1:sn-1]
    vρ1=vρ[1:sn-1]
    vθ2=vθ[sn:n]
    vρ2=vρ[sn:n]
    data=load_data(irflag,flag)
    cor1=cor(vρ1,vθ1)
    cor1=round(cor1,digits=2)
    cor2=cor(vρ2,vθ2)
    cor2=round(cor2,digits=2)
    Xρ=hcat(ones(n),data.ρ)
    regρt=ols(Xρ,vρ)
    betasρt=regρt.betas
    fρt(x)=betasρt[1].+betasρt[2].*x
    ρt=[fρt(x) for x ∈ minimum(data.ρ):0.001:maximum(data.ρ)]
    xgrid=minimum(data.ρ):0.001:maximum(data.ρ)
    p1=plot(xgrid,ρt,color=:red,label="Regression Line",xlabel="\\rho : Data",ylabel="\\rho : Model")
    scatter!(p1,data.ρ,vρ,color=:blue, label=false)
    p2=scatter(vρ1,vθ1,label ="Corr.Coef. before 1990=$cor1",ylim=(0,1))
    scatter!(vρ2,vθ2, label= "Corr.Coef. after 1990=$cor2", color=:red,shape=:x)
    p=plot(p1,p2,layout=(2,1))
return p
end

"""
figA2(flag)

Creates figure 5 for a specific country

"""
function figA2(flag::Int64)
    p1=figsubA2(1,flag,0)
    p2=figsubA2(1,flag,4)
    p3=figsubA2(1,flag,1)
    p=plot(p1,p2,p3,layout=(1,3),size=(1000,1000))
    return p
end

"""
figsubA3(flag)

Creates figure A3 for a specific country

"""
function figsubA3(flag::Int64)
    if flag==0
        title="Canada"
    elseif flag==1
        title="US"
    elseif flag==2
        title="Australia"
    else
        title="UK"
    end
    (l,m,vregime)=vs(1,flag,1)
    data=load_data(1,flag)
    p=scatter(data.year,vregime,xlabel="year",ylabel="Regime",color=:red,shape=:x,title="$title",ylim=(0.9,3.1))
    return p
end

"""
figA3(x)

Creates figure A3

"""
function figA3(x)
 p1= figsubA3(2)
 p3= figsubA3(0)
 p2= figsubA3(3)
 p4= figsubA3(1)
 p=plot(p1,p2,p3,p4,layout=(2,2),size=(1000,1000))
 return p
end


"""
figsubD2(irflag,flag,model)

Creates subplots of figure D2

"""
function figsubD2(irflag,flag,model)
    if flag==0
        title="Canada"
    elseif flag==1
        title="US"
    elseif flag==2
        title="Australia"
    else
        title="UK"
    end
    (vθ,vρ,vregime)=vs(irflag,flag,model)
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
    p1=plot(xgrid,ρt,color=:lightskyblue,label="Data Trend",xlabel="year",ylabel="\\rho",title="$title")
    scatter!(p1,data.year,data.ρ,color=:lightskyblue,label=false)
    plot!(p1,xgrid,vρt,color=:red,style=:dash,label="Model Trend")
    scatter!(p1,data.year,vρ,color=:red,shape=:cross,label=false)
end

"""
figD4(x)

Creates figure D2

"""
function figD2(x)
    p1=figsubD2(1,2,1)
    p2=figsubD2(2,2,1)
    p3=figsubD2(3,2,1)
    p4=figsubD2(1,0,1)
    p5=figsubD2(2,0,1)
    p6=figsubD2(3,0,1)
    p7=figsubD2(1,3,1)
    p8=figsubD2(2,3,1)
    p9=figsubD2(3,3,1)
    p10=figsubD2(1,1,1)
    p11=figsubD2(2,1,1)
    p12=figsubD2(3,1,1)

    p=plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,layout=(4,3),size=(2000,2000))
    return p
end

"""
figA4(x)

Creates figure A4

"""
function figA4(x)
    dict=matread("JiangShaoCodeData//atm_cic.mat")
    can=dict["ca_atm"]
    uk=dict["uk_atm"]
    usa=dict["us_atm"]
    aus=dict["aus_atm"]
    yr_can=can[:,1]
    fc_can=can[:,2]
    yr_uk=uk[:,1]
    fc_uk=uk[:,2]
    yr_usa=usa[:,1]
    fc_usa=usa[:,2]
    yr_aus=aus[:,1]
    fc_aus=aus[:,2]

    p1=plot(yr_aus,fc_aus,legend=false,xlabel="year",ylabel="ATM Withdrawal/CIC",title="Australia")
    p2=plot(yr_can,fc_can,legend=false,xlabel="year",ylabel="ATM Withdrawal/CIC",title="Canada")
    p3=plot(yr_uk,fc_uk,legend=false,xlabel="year",ylabel="ATM Withdrawal/CIC",title="UK")
    p4=plot(yr_usa,fc_usa,legend=false,xlabel="year",ylabel="ATM Withdrawal/CIC",title="US")

    p=plot(p1,p2,p3,p4,layout=(2,2),size=(1000,1000))

    return p
end

"""
figA5(x)

Creates figure A5

"""

function figA5(x)
    dict=matread("JiangShaoCodeData//atm_cic.mat")
    cr=dict["cr"]
    yr=cr[:,1]
    crc=cr[:,2]
    p=plot(yr,crc,xlabel="year", ylabel="Cash Receipts/CIC", size(500,500),legend=false)
 
    return p
end

"""
figD1(x)

Creates figure D1

"""
function figD1(x)
    dictir=matread("JiangShaoCodeData//ir.mat")
    usa=dictir["ir_usa"]
    yearusa=usa[:,1]
    ir1usa=usa[:,2]
    idx1=Vector{Bool}(undef,length(ir1usa))
    for i in 1:length(ir1usa)
        if (isnan(ir1usa[i]))
            idx1[i]=0
        else
            idx1[i]=1
        end
    end
    ir1usa=ir1usa[idx1]
    y1usa=yearusa[idx1]

    ir2usa=usa[:,3]
    idx2=Vector{Bool}(undef,length(ir2usa))
    for i in 1:length(ir2usa)
        if (isnan(ir2usa[i]))
            idx2[i]=0
        else
            idx2[i]=1
        end
    end
    ir2usa=ir2usa[idx2].*100
    y2usa=yearusa[idx2]

    ir3usa=usa[:,4].*100
    idx3=Vector{Bool}(undef,length(ir3usa))
    for i in 1:length(ir3usa)
        if (isnan(ir3usa[i]))
            idx3[i]=0
        else
            idx3[i]=1
        end
    end
    ir3usa=ir3usa[idx3]
    y3usa=yearusa[idx3]

    ir4usa=usa[:,5].*100
    idx4=Vector{Bool}(undef,length(ir4usa))
    for i in 1:length(ir4usa)
        if (isnan(ir4usa[i]))
            idx4[i]=0
        else
            idx4[i]=1
        end
    end
    ir4usa=ir4usa[idx4]
    y4usa=yearusa[idx4]

    ir5usa=usa[:,6].*100
    idx5=Vector{Bool}(undef,length(ir5usa))
    for i in 1:length(ir5usa)
        if (isnan(ir5usa[i]))
            idx5[i]=0
        else
            idx5[i]=1
        end
    end
    ir5usa=ir5usa[idx5]
    y5usa=yearusa[idx5]

    p1= plot(y1usa,ir1usa,color=:black, ylabel="Percentage",xlabel="year",title="Implied N. Interest Rate: US", label="Safe Bond")
    plot!(p1,y2usa,ir2usa,color=:red, style=:dashdot, label="DLM of Actual Inflation")
    plot!(p1,y3usa,ir3usa,color=:blue, style=:dot, label="Term Structure Based")
    plot!(p1,y4usa,ir4usa,color=:green, style=:dash, label="OECD Forecast")
    plot!(p1,y5usa,ir5usa,color=:lightskyblue, style=:dashdotdot, label="Survey Based")

    aus=dictir["ir_aus"]
    yearaus=aus[:,1]
    ir1aus=aus[:,2]
    idx1=Vector{Bool}(undef,length(ir1aus))
    for i in 1:length(ir1aus)
        if (isnan(ir1aus[i]))
            idx1[i]=0
        else
            idx1[i]=1
        end
    end
    ir1aus=ir1aus[idx1]
    y1aus=yearaus[idx1]

    ir2aus=aus[:,3]
    idx2=Vector{Bool}(undef,length(ir2aus))
    for i in 1:length(ir2aus)
        if (isnan(ir2aus[i]))
            idx2[i]=0
        else
            idx2[i]=1
        end
    end
    ir2aus=ir2aus[idx2].*100
    y2aus=yearaus[idx2]

    ir3aus=aus[:,4].*100
    idx3=Vector{Bool}(undef,length(ir3aus))
    for i in 1:length(ir3aus)
        if (isnan(ir3aus[i]))
            idx3[i]=0
        else
            idx3[i]=1
        end
    end
    ir3aus=ir3aus[idx3]
    y3aus=yearaus[idx3]

    ir4aus=aus[:,5].*100
    idx4=Vector{Bool}(undef,length(ir4aus))
    for i in 1:length(ir4aus)
        if (isnan(ir4aus[i]))
            idx4[i]=0
        else
            idx4[i]=1
        end
    end
    ir4aus=ir4aus[idx4]
    y4aus=yearaus[idx4]

    ir5aus=aus[:,6].*100
    idx5=Vector{Bool}(undef,length(ir5aus))
    for i in 1:length(ir5aus)
        if (isnan(ir5aus[i]))
            idx5[i]=0
        else
            idx5[i]=1
        end
    end
    ir5aus=ir5aus[idx5]
    y5aus=yearaus[idx5]

    p2= plot(y1aus,ir1aus,color=:black, ylabel="Percentage",xlabel="year",title="Implied N. Interest Rate: Australia", legend=false)
    plot!(p2,y2aus,ir2aus,color=:red, style=:dashdot)
    plot!(p2,y3aus,ir3aus,color=:blue, style=:dot)
    plot!(p2,y4aus,ir4aus,color=:green, style=:dash)
    plot!(p2,y5aus,ir5aus,color=:lightskyblue, style=:dashdotdot)

    can=dictir["ir_can"]
    yearcan=can[:,1]
    ir1can=can[:,2]
    idx1=Vector{Bool}(undef,length(ir1can))
    for i in 1:length(ir1can)
        if (isnan(ir1can[i]))
            idx1[i]=0
        else
            idx1[i]=1
        end
    end
    ir1can=ir1can[idx1]
    y1can=yearcan[idx1]

    ir2can=can[:,3]
    idx2=Vector{Bool}(undef,length(ir2can))
    for i in 1:length(ir2can)
        if (isnan(ir2can[i]))
            idx2[i]=0
        else
            idx2[i]=1
        end
    end
    ir2can=ir2can[idx2].*100
    y2can=yearcan[idx2]

    ir3can=can[:,4].*100
    idx3=Vector{Bool}(undef,length(ir3can))
    for i in 1:length(ir3can)
        if (isnan(ir3can[i]))
            idx3[i]=0
        else
            idx3[i]=1
        end
    end
    ir3can=ir3can[idx3]
    y3can=yearcan[idx3]

    ir4can=can[:,5].*100
    idx4=Vector{Bool}(undef,length(ir4can))
    for i in 1:length(ir4can)
        if (isnan(ir4can[i]))
            idx4[i]=0
        else
            idx4[i]=1
        end
    end
    ir4can=ir4can[idx4]
    y4can=yearcan[idx4]


    p3= plot(y1can,ir1can,color=:black, ylabel="Percentage",xlabel="year",title="Implied N. Interest Rate: Canada", legend=false)
    plot!(p3,y2can,ir2can,color=:red, style=:dashdot)
    plot!(p3,y3can,ir3can,color=:blue, style=:dot)
    plot!(p3,y4can,ir4can,color=:green, style=:dash)

    uk=dictir["ir_uk"]
    yearuk=uk[:,1]
    ir1uk=uk[:,2]
    idx1=Vector{Bool}(undef,length(ir1uk))
    for i in 1:length(ir1uk)
        if (isnan(ir1uk[i]))
            idx1[i]=0
        else
            idx1[i]=1
        end
    end
    ir1uk=ir1uk[idx1]
    y1uk=yearuk[idx1]

    ir2uk=uk[:,3]
    idx2=Vector{Bool}(undef,length(ir2uk))
    for i in 1:length(ir2uk)
        if (isnan(ir2uk[i]))
            idx2[i]=0
        else
            idx2[i]=1
        end
    end
    ir2uk=ir2uk[idx2].*100
    y2uk=yearuk[idx2]

    ir3uk=uk[:,4].*100
    idx3=Vector{Bool}(undef,length(ir3uk))
    for i in 1:length(ir3uk)
        if (isnan(ir3uk[i]))
            idx3[i]=0
        else
            idx3[i]=1
        end
    end
    ir3uk=ir3uk[idx3]
    y3uk=yearuk[idx3]

    ir4uk=uk[:,5].*100
    idx4=Vector{Bool}(undef,length(ir4uk))
    for i in 1:length(ir4uk)
        if (isnan(ir4uk[i]))
            idx4[i]=0
        else
            idx4[i]=1
        end
    end
    ir4uk=ir4uk[idx4]
    y4uk=yearuk[idx4]


    p4= plot(y1uk,ir1uk,color=:black, ylabel="Percentage",xlabel="year",title="Implied N. Interest Rate: UK", legend=false)
    plot!(p4,y2uk,ir2uk,color=:red, style=:dashdot)
    plot!(p4,y3uk,ir3uk,color=:blue, style=:dot)
    plot!(p4,y4uk,ir4uk,color=:green, style=:dash)

    p=plot(p1,p2,p3,p4,layout=(2,2),size=(1000,1000))

    return p
end


"""
Calibrate(irflag,flag,model)

Gives the value for parameters given a interest rate source, a country, and a model

"""
function Calibrate(irflag::Int64, flag::Int64, model::Int64)
    opt=[[0.0,0.0,0.0],0.0]
    if model==1
        if flag==3
            x0=[0.9 2.0 10.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        
        elseif (flag==1 && irflag==3)
            x0=[0.9,1.1,3.0]
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
        opt=Optim.optimize(x-> eqn_tfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        year2005=findall(x->x==2005,data.year)[1]
        year2014=findall(x->x==2014,data.year)[1]
        year1990=findall(x->x==1990,data.year)[1]
        θ=mean(vθ[year2005:year2014])
        cor_ρ=cor(vρ,data.ρ)[1]
        cor_θbef=cor(vρ[1:year1990-1],vθ[1:year1990-1])
        cor_θaft=cor(vρ[year1990:end],vθ[year1990:end])



        return (α,s,Xstar,rmsd,cor_ρ,cor_θbef,cor_θaft,qstar,θ)
    
    elseif model==0
        if flag==2
            x0=[0.5 1.0 15.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.9,1.0,10.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.0,0.9999999,0.0];
        ub=[0.999999,1.000001,Inf];
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vregime=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt=Optim.optimize(x-> eqn_LWfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        year2005=findall(x->x==2005,data.year)[1]
        year2014=findall(x->x==2014,data.year)[1]
        year1990=findall(x->x==1990,data.year)[1]
        θ=mean(vθ[year2005:year2014])
        cor_ρ=cor(vρ,data.ρ)[1]
        cor_θbef=cor(vρ[1:year1990-1],vθ[1:year1990-1])
        cor_θaft=cor(vρ[year1990:end],vθ[year1990:end])

        return (α,s,Xstar,rmsd,cor_ρ,cor_θbef,cor_θaft,qstar,θ)    

    elseif model==4
        if flag==3
            x0=[0.5 1.5 100.0]; #[eta s Xstar] initial guess
            x0=Vector{Float64}(vec(x0))
        else
            x0=[0.8,2.0,100.0];
            x0=Vector{Float64}(vec(x0))
        end
        lb=[0.001,0.001,0.0];
        ub=[1.0,Inf,Inf]
        R=[]
        r=[]
        data=load_data(irflag,flag)
        n=length(data.year)
        vz=zeros(n,1);
        vregime=zeros(n,1)
        vGDP=zeros(n,1);
        vzT=zeros(n,1)
        vρ=zeros(n,1)
        vθ=zeros(n,1)
        yi3=zeros(n,1)
        qi3=zeros(n,1)
        ystar=1.0
        vδprime=vδ(irflag,flag)
        opt=Optim.optimize(x-> eqn_noconnfit(x,irflag,flag),lb,ub,x0,Fminbox())
        x=opt.minimizer
        ess=opt.minimum
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
        year2005=findall(x->x==2005,data.year)[1]
        year2014=findall(x->x==2014,data.year)[1]
        year1990=findall(x->x==1990,data.year)[1]
        θ=mean(vθ[year2005:year2014])
        cor_ρ=cor(vρ,data.ρ)[1]
        cor_θbef=cor(vρ[1:year1990-1],vθ[1:year1990-1])
        cor_θaft=cor(vρ[year1990:end],vθ[year1990:end])

        return (α,s,Xstar,rmsd,cor_ρ,cor_θbef,cor_θaft,qstar,θ)
    end  
end

"""
Table1(x)

Store all results from Table 1 in a DataFrame

"""
function Table1(x)
  df=DataFrame()
  df.Country=["Canada","Canada","Canada","US","US","US","Australia","Australia","Australia","UK","UK","UK"]
  df.Model=["Standard", "NCF", "Full","Standard", "NCF", "Full","Standard", "NCF", "Full","Standard", "NCF", "Full"]
  α=zeros(12)
  s=zeros(12)
  Xstar=zeros(12)
  rmsd=zeros(12)
  cor_rho=zeros(12)
  cor_theta1=zeros(12)
  cor_theta2=zeros(12)
  qstar=zeros(12)
  θ=zeros(12)
    for mp in 1:3
        if mp==1
            model=0
        elseif mp==2
            model=4
        else
            model=1
        end
        for country in 0:3
            (α[country*3+mp],s[country*3+mp],Xstar[country*3+mp],rmsd[country*3+mp],cor_rho[country*3+mp],cor_theta1[country*3+mp],cor_theta2[country*3+mp],qstar[country*3+mp],θ[country*3+mp])=Calibrate(1,country,model)
        end 
    end
    df.α=α
    df.s=s
    df.xstar=Xstar
    df.NRMSE=rmsd
    df.cor_ρ=cor_rho
    df.cor_ρθ_before_1990=cor_theta1
    df.cor_ρθ_after_1990=cor_theta2
    df.Implied_qstar=qstar
    df.θ=θ
    return df
end

export  load_data,
        Calibrate,
        Table1,
        eqn_LWfit,
        eqn_Regime201610,
        eqn_noconnfit,
        eqn_tfit,
        fig5,
        figA2,
        figA3,
        figA4,
        figD1,
        figD4,
        figsub5,
        figsubA2,
        figsubA3,
        figsubD2,
        vs,
        vδ

end #module