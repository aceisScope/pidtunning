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
plot(cputime,current_t,'r',cputime,detector_t,'b');
