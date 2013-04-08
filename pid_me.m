c = csv_read('cooling_P100_I100_D1000.csv');
data = c{1,1};
setpoints = [];
errors = [];
setpoints = data(:,2);
errors = data(:,3);