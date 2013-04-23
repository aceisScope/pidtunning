clear all;
clc;

fid=fopen('warmup_20130419_13.58.50.csv');
c = textscan(fid,'%s %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid);

data = c{1,2};
cputime = data(:,2);

current_t = data (:,3);
detector_t = data (:,4);

figure(1);
title('Sine Function');
plot(cputime,current_t,'r',cputime,detector_t,'b');
xlabel('time');
ylabel('temperature');
axis([cputime(1),cputime(end),0,detector_t(end)]);

