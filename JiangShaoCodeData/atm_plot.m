clear;
figure(1);
load atm_cic.mat
subplot(2,2,1)
yr=aus_atm(:,1); % before 2005, data from Janet
fc=aus_atm(:,2); % after 2005, CPSS red book
plot(yr,fc);
ylabel('ATM Withdrawal/CIC')
xlabel('year')
title('Australia');

subplot(2,2,2)
yr=ca_atm(:,1);
fc=ca_atm(:,2);
plot(yr,fc);
ylabel('ATM Withdrawal/CIC')
xlabel('year')
title('Canada');

subplot(2,2,3)
yr=uk_atm(:,1);
fc=uk_atm(:,2);
plot(yr,fc);
ylabel('ATM Withdrawal/CIC')
xlabel('year')
title('UK');

subplot(2,2,4)
yr=us_atm(:,1);
fc=us_atm(:,2);
plot(yr,fc);
ylabel('ATM Withdrawal/CIC')
xlabel('year')
title('US');

figure(2)
yr=cr(:,1);
crc=cr(:,2);
plot(yr,crc);
ylabel('Cash Receipts/CIC')
xlabel('year')


