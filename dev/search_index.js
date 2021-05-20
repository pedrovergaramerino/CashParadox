var documenterSearchIndex = {"docs":
[{"location":"#Our-Replication-of-The-Cash-Paradox-(Jiang-and-Shao,-2019)","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"This replication study was part of our evaluation for the course Numerical Methods at SciencesPo Paris in Spring 2021","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"The functions used to replicate this paper are:","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"Modules = [CashParadox]","category":"page"},{"location":"#CashParadox.load_data","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.load_data","text":"load_data(irflag,flag)\n\nCreates a struct with the data for a given country and interest rate specification.\n\n\n\n\n\n","category":"type"},{"location":"#CashParadox.Calibrate-Tuple{Int64, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.Calibrate","text":"Calibrate(irflag,flag,model)\n\nGives the value for parameters given a interest rate source, a country, and a model.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.Table1-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.Table1","text":"Table1()\n\nStore all results from Table 1 in a DataFrame\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.Table2-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.Table2","text":"Table2(x)\n\nStore all results from Table 1 in a DataFrame\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.eqn_LWfit-Tuple{Vector{Float64}, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.eqn_LWfit","text":"eqn_LWfit\n\nMatching the time series of CIC_GDP ratio in the data with the simulated  series in the Lagos-Wright model.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.eqn_Regime201610-NTuple{6, Any}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.eqn_Regime201610","text":"eqn_Regime201610\n\nExtensive margin: This equation is necessary to solve for q2 and y in regime 2.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.eqn_noconnfit-Tuple{Any, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.eqn_noconnfit","text":"eqn_noconnfit\n\nMatching the time series of CIC_GDP ratio in the data with the simulated  series in the NCF model.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.eqn_tfit-Tuple{Any, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.eqn_tfit","text":"eqn_tfit\n\nMatching the time series of CIC_GDP ratio in the data with the simulated  series in the full model\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.fig5-Tuple{Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.fig5","text":"fig5(flag)\n\nCreates figure 5 for a specific country.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figA2-Tuple{Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figA2","text":"figA2(flag)\n\nCreates figure A2 for a specific country.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figA3-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figA3","text":"figA3()\n\nCreates figure A3 (all countries)\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figA4-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figA4","text":"figA4()\n\nCreates figure A4\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figD1-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figD1","text":"figD1(x)\n\nCreates figure D1\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figD2-Tuple{}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figD2","text":"figD2()\n\nReplicates figure D2.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figsub5-Tuple{Int64, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figsub5","text":"figsub5(irflag,flag,model)\n\nCreates the subplots for figure 5 for a given interest rate specification, model, and country.\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figsubA2-Tuple{Int64, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figsubA2","text":"figsubA2(irflag,flag,model)\n\nCreates the subplots for figure A2 for a given interest rate specification, model, and country\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figsubA3-Tuple{Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figsubA3","text":"figsubA3(flag)\n\nCreates figure A3 for a specific country\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.figsubD2-Tuple{Any, Any, Any}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.figsubD2","text":"figsubD2(irflag,flag,model)\n\nCreates subplots of figure D2 for a given country, interest rate specification and model\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.vs-Tuple{Int64, Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.vs","text":"vs(irflag,flag,model)\n\nReturns simulated θ and simulated ρ for a specified model\n\n\n\n\n\n","category":"method"},{"location":"#CashParadox.vδ-Tuple{Int64, Int64}","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"CashParadox.vδ","text":"vδ(irflag,flag)\n\nComputes vδ given a country and interest rate specification.\n\n\n\n\n\n","category":"method"},{"location":"#Replication-of-model-predictions:-Figures-5a-5d","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of model predictions: Figures 5a-5d","text":"","category":"section"},{"location":"#Figure-5a:-AUSTRALIA","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure 5a: AUSTRALIA","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Fig5(2)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: fig5aus)","category":"page"},{"location":"#Figure-5b:-CANADA","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure 5b: CANADA","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Fig5(0)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: fig5can)","category":"page"},{"location":"#Figure-5c:-UK","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure 5c:  UK","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Fig5(3)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: fig5uk)","category":"page"},{"location":"#Figure-5d:-US","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure 5d: US","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Fig5(1)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: fig5us)","category":"page"},{"location":"#Replication-of-model-predictions:-Figures-A2a-A2d","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of model predictions: Figures A2a-A2d","text":"","category":"section"},{"location":"#Figure-A2a:-AUSTRALIA","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure A2a: AUSTRALIA","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA2(2)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA2aus)","category":"page"},{"location":"#Figure-A2b:-CANADA","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure A2b: CANADA","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA2(0)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA2can)","category":"page"},{"location":"#Figure-A2c:-UK","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure A2c:  UK","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA2(3)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA2uk)","category":"page"},{"location":"#Figure-A2d:-US","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Figure A2d: US","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA2(1)","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA2us)","category":"page"},{"location":"#Replication-of-regime-changes-over-time:-Figure-A3","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of regime changes over time: Figure A3","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA3()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA3)","category":"page"},{"location":"#Replication-of-The-Value-of-ATM-withdrawals-over-CIC:-Figure-A4","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of The Value of ATM withdrawals over CIC: Figure A4","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA4()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA4)","category":"page"},{"location":"#Replication-of-Cash-receipts-from-circulation-in-the-Federal-Reserve-Banks:-Figure-A5","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of Cash receipts from circulation in the Federal Reserve Banks: Figure A5","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigA5()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figA5)","category":"page"},{"location":"#Replication-of-Different-measures-of-nominal-interest-rates:-Figure-D1","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of Different measures of nominal interest rates: Figure D1","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigD1()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figD1)","category":"page"},{"location":"#Replication-of-CIC/GDP-with-different-interest-rate-specifications:-Figure-D2","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of CIC/GDP with different interest rate specifications: Figure D2","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create this figure, type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> FigD2()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figD2)","category":"page"},{"location":"#Replication-of-Table-1-and-Table-2:-Calibration-results-and-Cash-shares-relative-to-credit","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of Table 1 and Table 2: Calibration results and Cash shares relative to credit","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"Most of results in this table are very close to those from Jiang & Shao (2019). However, we encounter problems replicating the result for the NCF model with UK data.","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create a dataframe containing this results type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Table1()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figG1)","category":"page"},{"location":"#Replication-of-Table-D.1:-Parameter-values-and-Table-D.2:-Model-performance-comparison","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Replication of Table D.1: Parameter values and Table D.2: Model performance comparison","text":"","category":"section"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"In order to create a dataframe containing this results type","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"julia> using Cash Paradox\n\njulia> Table2()","category":"page"},{"location":"","page":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","title":"Our Replication of The Cash Paradox (Jiang & Shao, 2019)","text":"(Image: figD1D2)","category":"page"}]
}
