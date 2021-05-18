These Matlab codes are used to replicate the results in the paper "The Cash Paradox" by Janet Jiang and Enchuan Shao.

Data file:
data.mat: contains prime corporate bond yields and cash/GDP ratios
ir.mat: contains different nominal interest rates computed by different measures of expected inflation
atm_cic.mat: contains the time series of ATM withdrawals over CIC for different countries, and cash receipts from the Fed.

Program file:
calibrate_tfit.m: main program to calibrate and simulate different models, graphical results will be generated. The main program will call the following subroutines correspondingly:
eqn_Regime201610.m, eqn_tfit.m, eqn_tfit2.m, eqn_LWfit.m, eqn_noconnfit.m, eqn_noconnfit1.m, and eqn_noconnfit2.m.
1. To generate Figures 5(a) to 5(d), A.2(a) to A.2(d), set fig=1, then select country counter, flag=0, 1, 2, 3 one by one. 
2. To generate Figures A.3 and D.2, choose fig=2.  

plotid.m: plots the data on nominal interest rates and credit expansion rates (Figures A.1 and D.1).

atm_plot.m: plots Figures A.4 and A.5

Expected computation time: 1 minute.

