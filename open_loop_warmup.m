clear all;
clc;
clf;

fid=fopen('warmup_20130419_13.58.50.csv');
c1 = textscan(fid,'%s %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid);

data = c1{1,2};

cputime = data(:,2);
starttime = cputime (1);
cputime = cputime - starttime;

current_t = data (:,3); %cold finger
do = current_t(1);
current_t = current_t - do;

detector_t = data (:,4); %detector
tec_current2 = data (:,5);

figure(1);
title('open loop');
plot(cputime,current_t,'r');
xlabel('time(sec)');
ylabel('temperature(celcius)');
axis([cputime(1),cputime(end),current_t(1),current_t(end)]);

%two-point method
t2 = cputime(1000);
y2 = current_t(1000);
t1 = cputime (160);
y1 = current_t(160);
y_stable = current_t(end);
% K = y_stable - current_t(1);
K = y_stable;
T = (t2 - t1)/ (log(1-y1/y_stable) - log(1-y2/y_stable));
tau = (t2*log(1-y1/y_stable) - t1*log(1-y2/y_stable))/(log(1-y1/y_stable) - log(1-y2/y_stable)); 


%estimate
K = y_stable - current_t(1);
tau = 1.063;
T = cputime(204) - tau; %where y(t)= y(inf)*63.2%

check = [];
for i = 1:length(cputime)
    check = [check, K*(1-exp((-(cputime(i)-tau))/T))];
end
% check = check - do;
hold on;
plot(cputime,check,'g');


%fitcurve
start_point = [70 210 1]';
[estimates,mode1]=fitcurve(cputime,current_t,start_point);

K = estimates(1);
T = estimates(2);
if estimates(3) >= 0.5*tau && estimates(3) <= 1.5*tau
    tau = estimates(3);
end


check = [];
for i = 1:length(cputime)
    check = [check, K*(1-exp((-(cputime(i)-tau))/T))];
end
% check = check - do;
hold on;
plot(cputime,check,'black');

%step of tf
sys = tf([K],[T 1],'ioDelay',tau);
step(sys); 

legend('current temperature','engineering estimate','fit curve');
