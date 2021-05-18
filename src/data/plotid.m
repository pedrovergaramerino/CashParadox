clear;
close;
set(gcf,'DefaultLineLineWidth',1.0);
load data;

f1=figure(1);
% Australia
% load ir_aus.txt;
% year=ir_aus(:,1);
% ir=ir_aus(:,2);
year=aus(:,1);
ir=aus(:,2);
ddata=[0,0.520354595,0.537666422,0.57032816,0.594028535,0.622316221,...
    0.65072652,0.687205136,0.694317256,0.73605656,0.760156768,0.792611398,...
    0.841961623,0.878037024,0.905946326,0.90592637,0.897747613,0.904361411,...
    0.90064194]';
yt=[1960,1994:2011]';
% Use the data to find the fitted curve of delta
[dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
subplot(4,2,1)
plot(yt,ddata,'o');
hold on
p1=plot(dt,'predobs');
ylabel('\delta')
xlabel('year')
xlim([1960 2015]);
hL=legend(p1,'fitted curve','95% Confidence Bounds','Orientation','horizontal');
title('Canada');
newPosition = [0.4 0.4 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

hold off
title('Credit Expansion: Australia')
subplot(4,2,2)
[ir_t,gof1]=fit(year,ir,'poly2','Normalize','on','Robust','LAR');
plot(year,ir,'o');
hold on
p1=plot(ir_t,'predobs');
ylabel('N. Int. Rate')
xlabel('year')
xlim([1960 2015]);
hold off
title('N. Interest Rate: Australia');
b=gca;
legend(b,'off') 

% Canada
% clear;
% load ir_can.txt;
% year=ir_can(:,1);
% ir=ir_can(:,2);
year=can(:,1);
ir=can(:,2);
ddata=[0,0.49,0.782,0.822,0.856]'; % Data on access to credit
yt=[1960,1980,1999,2005,2012]';
[dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
subplot(4,2,3)
plot(yt,ddata,'o');
hold on
plot(dt,'predobs');
ylabel('\delta')
xlabel('year')
xlim([1960 2015]);
hold off
b=gca;
legend(b,'off') 
title('Credit Expansion: Canada')
subplot(4,2,4)
[ir_t,gof1]=fit(year,ir,'poly2','Normalize','on','Robust','LAR');
plot(year,ir,'o');
hold on
plot(ir_t,'predobs');
ylabel('N. Int. Rate')
xlabel('year')
xlim([1960 2015]);
hold off
title('N. Interest Rate: Canada');
b=gca;
legend(b,'off') 

% UK
% clear;
% load ir_uk.txt;
% year=ir_uk(:,1);
% ir=ir_uk(:,2);
year=uk(:,1);
ir=uk(:,2);
ddata=[0,0.5,0.53,0.65,0.62,0.64,0.62,0.61, 0.6]';
yt=[1960,1999, 2000, 2005, 2009, 2010, 2011, 2012, 2013]';
[dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
subplot(4,2,5)
plot(yt,ddata,'o');
hold on
plot(dt,'predobs');
ylabel('\delta')
xlabel('year')
xlim([1960 2015]);
hold off
b=gca;
legend(b,'off') 
title('Credit Expansion: UK')
subplot(4,2,6)
[ir_t,gof1]=fit(year,ir,'poly2','Normalize','on','Robust','LAR');
plot(year,ir,'o');
hold on
plot(ir_t,'predobs');
ylabel('N. Int. Rate')
xlabel('year')
xlim([1960 2015]);
hold off
b=gca;
legend(b,'off') 
title('N. Interest Rate: UK');

% US
% clear;
% load ir_usa.txt;
% year=ir_usa(:,1);
% ir=ir_usa(:,2);
year=usa(:,1);
ir=usa(:,2);
ddata=[0,16.3,38.3,43.0,55.8,62.2,66.5,68,72.6,71.5]'/100;
yt=[1960,1970,1977,1983,1989,1992,1995,1998,2001,2004]';
[dt,gof]=fit(yt,ddata,'poly2','Normalize','on','Robust','LAR');
subplot(4,2,7)
plot(yt,ddata,'o');
hold on
plot(dt,'predobs');
ylabel('\delta')
xlabel('year')
xlim([1960 2015]);
hold off
b=gca;
legend(b,'off') 
title('Credit Expansion: US')
subplot(4,2,8)
[ir_t,gof1]=fit(year,ir,'poly2','Normalize','on','Robust','LAR');
plot(year,ir,'o');
hold on
plot(ir_t,'predobs');
ylabel('N. Int. Rate')
xlabel('year')
xlim([1960 2015]);
hold off
b=gca;
legend(b,'off') 
title('N. Interest Rate: US');

f2=figure(2);
clear;
set(gcf,'DefaultLineLineWidth',2.0);

load ir;
% US
subplot(2,2,1)
year=ir_usa(:,1);
ir1=ir_usa(:,2); % Corporate bond/commercial paper/bankers acceptance
idx=isnan(ir1);
y1=year(~idx);
ir1=ir1(~idx);
ir2=ir_usa(:,3);    % DLM of actual inflation
idx=isnan(ir2);
y2=year(~idx);
ir2=ir2(~idx)*100;
ir3=ir_usa(:,4);    % Term structure based
idx=isnan(ir3);
y3=year(~idx);
ir3=ir3(~idx)*100;
ir4=ir_usa(:,5);    % OECD Forecast
idx=isnan(ir4);
y4=year(~idx);
ir4=ir4(~idx)*100;  
ir5=ir_usa(:,6);    % Survey based
idx=isnan(ir5);
y5=year(~idx);
ir5=ir5(~idx)*100;  
p2=plot(y1,ir1,'-k',y2,ir2,'-.r',y3,ir3,':b',y4,ir4,'--g',y5,ir5,'oc');
title('Implied N. Interest Rate: US')
ylabel('Percentage')
xlabel('year')
hL=legend(p2,'Safe Bond','DLM of Actual Inflation','Term Structure Based','OECD Forecast','Survey Based','Orientation','horizontal');
newPosition = [0.4 0.4 0.2 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
% Australia
subplot(2,2,2)
year=ir_aus(:,1);
ir1=ir_aus(:,2); % Corporate bond/commercial paper/bankers acceptance
idx=isnan(ir1);
y1=year(~idx);
ir1=ir1(~idx);
ir2=ir_aus(:,3);    % DLM of actual inflation
idx=isnan(ir2);
y2=year(~idx);
ir2=ir2(~idx)*100;
ir3=ir_aus(:,4);    % Term structure based
idx=isnan(ir3);
y3=year(~idx);
ir3=ir3(~idx)*100;
ir4=ir_aus(:,5);    % OECD Forecast
idx=isnan(ir4);
y4=year(~idx);
ir4=ir4(~idx)*100;  
ir5=ir_aus(:,6);    % Survey based
idx=isnan(ir5);
y5=year(~idx);
ir5=ir5(~idx)*100;  
plot(y1,ir1,'-k',y2,ir2,'-.r',y3,ir3,':b',y4,ir4,'--g',y5,ir5,'oc');
title('Implied N. Interest Rate: Australia')
ylabel('Percentage')
xlabel('year')

% Canada
subplot(2,2,3)
year=ir_can(:,1);
ir1=ir_can(:,2); % Corporate bond/commercial paper/bankers acceptance
idx=isnan(ir1);
y1=year(~idx);
ir1=ir1(~idx);
ir2=ir_can(:,3);    % DLM of actual inflation
idx=isnan(ir2);
y2=year(~idx);
ir2=ir2(~idx)*100;
ir3=ir_can(:,4);    % Term structure based
idx=isnan(ir3);
y3=year(~idx);
ir3=ir3(~idx)*100;
ir4=ir_can(:,5);    % OECD Forecast
idx=isnan(ir4);
y4=year(~idx);
ir4=ir4(~idx)*100;  
plot(y1,ir1,'-k',y2,ir2,'-.r',y3,ir3,':b',y4,ir4,'--g');
title('Implied N. Interest Rate: Canada')
ylabel('Percentage')
xlabel('year')

% UK
subplot(2,2,4)
year=ir_uk(:,1);
ir1=ir_uk(:,2); % Corporate bond/commercial paper/bankers acceptance
idx=isnan(ir1);
y1=year(~idx);
ir1=ir1(~idx);
ir2=ir_uk(:,3);    % DLM of actual inflation
idx=isnan(ir2);
y2=year(~idx);
ir2=ir2(~idx)*100;
ir3=ir_uk(:,4);    % Term structure based
idx=isnan(ir3);
y3=year(~idx);
ir3=ir3(~idx)*100;
ir4=ir_uk(:,5);    % OECD Forecast
idx=isnan(ir4);
y4=year(~idx);
ir4=ir4(~idx)*100;  
plot(y1,ir1,'-k',y2,ir2,'-.r',y3,ir3,':b',y4,ir4,'--g');
title('Implied N. Interest Rate: UK')
ylabel('Percentage')
xlabel('year')


