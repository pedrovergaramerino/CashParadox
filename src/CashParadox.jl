module CashParadox
using MAT

"""
   load_data(fig,irflag,flag,)

Load data depending on the figure and country wanted
"""
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
            ddata=[0,0.49,0.782,0.822,0.856]' #Data on access to credit
            yt=[1960,1980,1999,2005,2012]'
            ρ=can[:,3]
        elseif flag==1 #USA
            year=usa[:,1]
            ir=usa[:,2]
            ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100
            yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]'
            ρ=usa[:,3]
        elseif flag==2 #AUSTRALIA
            year=aus[:,1]
            ir=aus[:,2]
            ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,0.90064194]'
            yt= vcat([1965], collect(1994:2011))
            ρ=aus[:,3]
        else #UK 
            uk=uk[5:end,:]
            year=uk[:,1]
            ir=uk[:,2]
            ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]'
            yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]'
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
            ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100
            yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]'
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
            ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,0.90064194]';
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
            ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]'
            yt=[1965,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]'
        end
    end
    data=Dict(:year => year, :ir => ir, :ρ => ρ, :ddata => ddata, :yt=> yt)
    return data
end


"""
   eqn_LWfit(x)

Matching the time series of CIC_GDP ratio in the data with the simulated series in the Lagos-Wright model. 
"""
function eqn_LWfit(x::Vector{Float64})
η=x[1]
s=x[2]
Xstar=x[3]
α=η
n=length(year)
end

"""
    domath(x::Number)

Return `x + 5`.
"""
domath(x::Number) = x + 5

end
