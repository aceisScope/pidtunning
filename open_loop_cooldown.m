clear all;
clc;
clf;

fid=fopen('cooldown_20130419_13.24.44.csv');
c1 = textscan(fid,'%s %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid);

data = c1{1,2};

cputime = data(:,2);
starttime = cputime (1);
cputime = cputime - starttime;

current_t = data (:,3); %cold finger
do = current_t(1);
current_t = abs(current_t - do);

detector_t = data (:,4); %detector
tec_current2 = data (:,5);

figure(1);
title('open loop');
plot(cputime,current_t,'r');
xlabel('time(sec)');
ylabel('temperature(celcius)');

y_stable = current_t(end);
K = abs(y_stable - current_t(1));
tau = 4.184; % 1366403091.216000 - 1366403088.076000
T = cputime(156) - tau; %where y(t)= y(inf)*63.2%

check = [];
for i = 1:length(cputime)
    check = [check, K*(1- exp(-(cputime(i)-tau)/T))];    
end
hold on;
plot(cputime,check,'g');

%fitcurve
start_point = [80 170 3]';
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
title('cooldown 20130419 13.24.44');

legend('current temperature','engineering estimate','fit curve');



fid=fopen('cooldown_20130419_15.45.06.csv');
c1 = textscan(fid,'%s %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid);

data = c1{1,2};

cputime = data(:,2);
starttime = cputime (1);
cputime = cputime - starttime;

current_t = data (:,3); %cold finger
do = current_t(1);
current_t = abs(current_t - do);

detector_t = data (:,4); %detector
tec_current2 = data (:,5);

figure(2);
title('open loop');
plot(cputime,current_t,'r');
xlabel('time(sec)');
ylabel('temperature(celcius)');

y_stable = current_t(end);
K = abs(y_stable - current_t(1));
tau = 3.198; % 1366411513.907000 - 1366411510.709000
T = cputime(166) - tau; %where y(t)= y(inf)*63.2%

check = [];
for i = 1:length(cputime)
    check = [check, K*(1- exp(-(cputime(i)-tau)/T))];    
end
hold on;
plot(cputime,check,'g');

%fitcurve
start_point = [80 170 3]';
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
title('cooldown 20130419 15.45.06');

legend('current temperature','engineering estimate','fit curve');