function csv_data = csv_result(filename)
fid=fopen(filename);
% fid=fopen('cooling_P100_I100_D1000.csv');
csv_data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','HeaderLines',1,'Delimiter',',','CollectOutput',1);
fclose(fid);
 


