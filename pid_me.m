clear all;
clc;

%discrete PID simulation
c = csv_read('cooling_P100_I100_D1000.csv'); %reading data from csv into cell
data = c{1,1};
inits = data(1,:);
Errors = data(:,3);
CurrentMeasurments = data(:,10);

%external parametres
timestamp = 0;
setpoint = inits(2); 
temperature = inits(3); 
error = inits(4);
sum_error = inits(5);
proportionalFactor = 0;
integralFactor = 0;
derivativeFactor = 0; 
control = inits(9);
lastCurrentMeasurment = CurrentMeasurments(1);

%PID practice...
Ki = 100;
Kp = 100;
Kd = 1000; 

%Sampling time
ts = 0.01;

controls(1) = control;
time(1) = ts*1;

%control range [0,60100] (100uA)
br_MaxPidControlValue = 60100.0;
	
for k=2:1:847
time(k) = k*ts;

error_last = error;
error = Errors(k);

sum_error = sum_error + error;

%Clip "sum" to +/- 'br_MaxPidControlValue'
if ( br_MaxPidControlValue < sum_error)
    sum_error = br_MaxPidControlValue;
end;
if ( -br_MaxPidControlValue > sum_error)
    sum_error = -br_MaxPidControlValue;
end;

control_last = control;

%proportional 
proportionalFactor = Kp*error;
control = proportionalFactor;

%integral 
integralFactor = Ki*sum_error;
control = control + integralFactor;

%derivative 
derivativeFactor = Kd*(error-error_last);
control = control + derivativeFactor;

% Limit increases in current to 1/16 of maximum.
% Limit absolute current to maximum for the control (60100 * 100uA = 6.01 Amps).
% Limit minimum current to 0.
%Clip "control" to 0 and 'br_MaxPidControlValue'
stepsize = 	br_MaxPidControlValue*0.001;
if (stepsize < (control - control_last)) 
    control = control + stepsize;
end;
if ( br_MaxPidControlValue < control)
    control = br_MaxPidControlValue;
end;
if ( 0 > control)
    control = 0;
end;

controls(k) = control;
	
end;

figure(1);
plot(time,CurrentMeasurments,'r',time,controls,'b');
	

