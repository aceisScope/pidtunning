clear all;
clc;
clf;

fid=fopen('warmup_20130419_13.58.50.csv');
c = textscan(fid,'%s %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid);

data = c{1,2};
cputime = data(:,2);
starttime = cputime (1);
cputime = cputime - starttime;
current_t = data (:,3); %cold finger
detector_t = data (:,4); %dector
tec_current = data (:,5);

figure(1);
title('open loop');
plot(cputime,current_t,'r',cputime,detector_t,'b');
xlabel('time(sec)');
ylabel('temperature(celcius)');
axis([cputime(1),cputime(end),current_t(1),current_t(end)]);


t2 = cputime(1000);
y2 = current_t(1000);
t1 = cputime (160);
y1 = current_t(160);
y_stable = current_t(end);
K = (y_stable - current_t(1))/tec_current(1);
T = (t2 - t1)/ (log(1-y1/y_stable) - log(1-y2/y_stable));
tau = (t2*log(1-y1/y_stable) - t1*log(1-y2/y_stable))/(log(1-y1/y_stable) - log(1-y2/y_stable)); 


%estimate
K = y_stable - current_t(1);
tau = 1.063;
T = cputime(204); %where y(t)= y(inf)*63.2%

%fitcurve
start_point = [K T tau current_t(1)];
[estimates,mode1]=fitcurve(cputime,current_t,start_point);

check = [];
for i = 1:length(cputime)
    check = [check, K*(1-exp((-(cputime(i)-tau))/T)) - 44.10];
end
hold on;
plot(cputime,check,'g');


