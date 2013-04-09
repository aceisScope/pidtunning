%discrete PID simulation
c = csv_read('cooling_P100_I100_D1000.csv'); %reading data from csv into cell
data = c{1,1};
inits = data(0,:);
Errors = data(:,3);
CurrentMeasurments = data(:,0);

%external parametres
timestamp = 0;
setpoint = inits(1); 
temperature = inits(2); 
error = inits(3);
sum = inits(4);
proportionalFactor = 0;
integralFactor = 0;
derivativeFactor = 0; 
control = inits(8);
lastCurrentMeasurment = CurrentMeasurements(0);

%PID practice...
Ki = 100;
Kp = 100;
Kd = 1000; 

%Sampling time
ts = 0.01s;

%control range [0,60100] (100uA)
br_MaxPidControlValue = 60100.0f;
	
for k=1:1:2000
	
	time(k) = k*ts;

	error_last = error;
	error = Errors(k);
	
	sum += error;
	
	%Clip "sum" to +/- 'br_MaxPidControlValue'
	if ( br_MaxPidControlValue < sum)
		sum = br_MaxPidControlValue;
	end;
	if ( -br_MaxPidControlValue > sum)
		sum = -br_MaxPidControlValue;
	end;

	control_last = control;
	
	%proportional 
	proportionalFactor = Kp*error;
	control = proportionalFactor;
	
	%integral 
	integralFactor = ki*sum;
	control += integralFactor;
		
	%derivative 
	derivativeFactor = kd*(error-error_last);
	control += derivativeFactor;
	
	% Limit increases in current to 1/16 of maximum.
	% Limit absolute current to maximum for the control (60100 * 100uA = 6.01 Amps).
	% Limit minimum current to 0.
	%Clip "control" to 0 and 'br_MaxPidControlValue'
	stepsize = 	br_MaxPidControlValue*0.001;
	if (stepsize < (control - control_last)) 
		control += stepsize;
	end;
	if ( br_MaxPidControlValue < sum)
		control = br_MaxPidControlValue;
	end;
	if ( 0 > control)
		control = 0;
	end;
	
end;
	

