% Copyright (C) 2002: Juergen Hahn
% Department of Chemical Engineering
% Texas A&M University
% College Station, TX 77843-3122
% USA
% hahn@che.tamu.edu
%
% This code is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% This code is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details at 
% http://www.gnu.org/licenses/gpl.html


function varargout = Controller_Design(varargin)
% CONTROLLER_DESIGN Application M-file for Controller_Design.fig
%    FIG = CONTROLLER_DESIGN launch Controller_Design GUI.
%    CONTROLLER_DESIGN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 23-Sep-2002 10:41:56


global y t model_order func prop integ deriv G_c sys_fit sys Kc tau_c tau_i tau_d prop integ deriv K theta tau_1 tau_2 beta integrating G_c_non_ideal data_imported

if isempty(model_order) == 1
    model_order = 1;
else
end

if isempty(tau_c) == 1
    tau_c = 1.0;
else
end

if isempty(y) == 1
    % Gain
    K = 1;
    % dominant time constant
    tau_1 = 1;
    % other time constant
    tau_2 = 0.2;
    % Time delay
    theta = 0.5;
        
    sys_fit = tf([K],conv([tau_1 1],[tau_2 1]));
    sys_fit.inputdelay = theta;
    sys     = sys_fit;
    [y,t]   = step(sys_fit);
    
else
end

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    if model_order == 1
        set(handles.PI,'Value',1);
        set(handles.PID,'Value',0);
    else
        set(handles.PI,'Value',0);
        set(handles.PID,'Value',1);
    end
    
    if tau_c ~= 1.0
        set(handles.enter_tau_c,'String',num2str(tau_c));
    else
    end

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = PI_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys

set(handles.PID,'Value',0);

model_order = 1;


% --------------------------------------------------------------------
function varargout = PID_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys

set(handles.PI,'Value',0);

model_order = 2;


% --------------------------------------------------------------------
function varargout = Compute_Controller_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys Kc tau_c tau_i tau_d prop integ deriv K theta tau_1 tau_2 beta integrating G_c_non_ideal data_imported tau chsi

% controller parameter bounds

tau_c_min = 0.05;

% factor for determining tau_c (I would usually recommend fac = 1, but for some examples higher values might be useful)
% fac = 1 means that closed-loop time constant can only be at most 8 times smaller than dominant open-loop time constant

fac = 5.0;

% Pre-process data
i = 1;
theta_initial = t(i);
while y(i) == y(i+1),
    i = i + 1;
    theta_initial = t(i);
end

% Identify transfer function model from data

% This will depend upon if you want a PI or a PID controller
% func = 1; % PI from first order process plus time delay
% func = 2; % PID from second order process plus time delay
% func = 3; % pure I from pure time delay process
% func = 4; % PI from integrating process plus time delay
% func = 5; % PID from integrating process plus time delay plus first order process
% func = 6; % PID from double integrating process plus time delay
% func = 7; % PI from process with inverse response plus time delay
% func = 8; % PI from integrating process plus time delay plus inverse response

switch model_order
    
    % For a PI controller a first order transfer function needs to be chosen
    case {1}

        K_initial = max(y);
        tau_initial = t(length(t))/5;
        beta_initial = t(length(t))/5;
    
        x_initial = [theta_initial K_initial tau_initial];

        x = x_initial;
        
        disp('Computing model structure and parameters for PI controller.');
     
        [x_1, resnorm_1] = lsqnonlin(@first_order_fit_function,x,[0 min(y) 0],[t(length(t)) max(y) t(length(t))],[],y,t);
        if y(length(y)) >= y(length(y)-1)
            [x_2, resnorm_2] = lsqnonlin(@first_order_fit_function_2,[theta_initial x_1(2)/x_1(3)],[0 0],[t(length(t)) max([max(abs(y))/t(length(t))*3 max(abs(y))])],[],y,t);
        else
            [x_2, resnorm_2] = lsqnonlin(@first_order_fit_function_2,[theta_initial x_1(2)/x_1(3)],[0 -max([max(abs(y))/t(length(t))*3 max(abs(y))])],[t(length(t)) 0],[],y,t);
        end
        [x_3, resnorm_3] = lsqnonlin(@first_order_fit_function_3,[theta_initial x_1(2)],[theta_initial min(y)],[t(i+1) max(y)],[],y,t);
        [x_4, resnorm_4] = lsqnonlin(@first_order_fit_function_4,[x_1 beta_initial],[t(i) min(y) 0 0],[t(length(t)) max(y) t(length(t)) t(length(t))*20],[],y,t);
        if y(length(y)) >= y(length(y)-1)
            [x_5, resnorm_5] = lsqnonlin(@first_order_fit_function_5,[theta_initial (y(length(y))-y(i+1))/(t(length(t))-t(i+1)) -y(i+1)/((y(length(y))-y(i+1))/(t(length(t))-t(i+1)))],[t(i) 0 0],[t(i+1) max(abs(y))/t(length(t)) t(length(t))*20],[],y,t);
        else
            [x_5, resnorm_5] = lsqnonlin(@first_order_fit_function_5,[theta_initial (y(length(y))-y(i+1))/(t(length(t))-t(i+1)) -y(i+1)/((y(length(y))-y(i+1))/(t(length(t))-t(i+1)))],[t(i) -max(abs(y))/t(length(t)) 0],[t(i+1) 0 t(length(t))*20],[],y,t);
        end
        %resnorm_5 = 1E16;
        % the last case seems to cause numerical problem and never converges unless you are exactly at the correct solution
        resnorm = min([resnorm_1, resnorm_2, resnorm_3, resnorm_4, resnorm_5]);
                      
        if resnorm_3 == resnorm
            % Pure time delay
            K       = x_3(2);
            tau_1   = 0;
            tau_2   = 0;
            theta   = x_3(1);
            beta    = 0;
            sys_fit = tf([K],[1]);
            sys_fit.inputdelay = theta;
            tau_c = max([theta, tau_c_min, t(length(t))/(40.0*fac)]);
            func    = 3;
            integrating = 0;
        else
        end
        if resnorm_2 == resnorm
            % Integrating plus time delay
            K       = x_2(2);
            tau_1   = 0;
            tau_2   = 0;
            theta   = x_2(1);
            beta    = 0;
            sys_fit = tf([K],[1 0]);
            sys_fit.inputdelay = theta;
            tau_c = max([theta, tau_c_min, t(length(t))/(10.0*fac)]);
            func    = 4;
            integrating = 1;
        else
        end
        if resnorm_5 == resnorm
            % Integrating plus time delay plus inverse response
            K       = x_5(2);
            tau_1   = 0;
            tau_2   = 0;
            theta   = x_5(1);
            beta    = x_5(3);
            sys_fit = tf(K*[-beta 1],[1 0]);
            sys_fit.inputdelay = theta;
            tau_c = max([theta, tau_c_min, beta+theta, t(length(t))/(40.0*fac)]);
            func    = 8;
            integrating = 1;
        else
        end
        if resnorm_1 == resnorm
            % First order transfer function plus time delay
            K       = x_1(2);
            tau_1   = x_1(3);
            tau_2   = 0;
            theta   = x_1(1);
            beta    = 0;
            sys_fit = tf([K],[tau_1 1]);
            sys_fit.inputdelay = theta;
            if theta * 8 < tau_1
                tau_c = max([theta, tau_1/(8*fac), tau_c_min, t(length(t))/(40.0*fac)]);
            else
                tau_c = max([theta, tau_c_min, t(length(t))/(40.0*fac)]);
            end
            func    = 1;
            integrating = 0;
        else
        end
        if resnorm_4 == resnorm
            % First order transfer function plus time delay plus inverse response
            K       = x_4(2);
            tau_1   = x_4(3);
            tau_2   = 0;
            theta   = x_4(1);
            beta    = x_4(4);
            sys_fit = tf(K*[-beta 1],[tau_1 1]);
            sys_fit.inputdelay = theta;
            if theta * 8 < tau_1
                tau_c = max([theta, tau_1/(8*fac), tau_c_min, beta+theta, t(length(t))/(40.0*fac)]);
            else
                tau_c = max([theta, tau_c_min, beta, t(length(t))/(40.0*fac)]);
            end
            func    = 7;
            integrating = 0;
        else
        end
                
        [y_fit,t_fit] = step(sys_fit,t);
        axes(handles.axes1);
        plot(t,y,'k-',t_fit,y_fit,'k:');
        legend('Data','First-order fit');
        xlabel('t');
        ylabel('y');
        title('y vs. t for original model and fitted model that is used for controller design');
        clear x_1 x_2 x_3 resnorm_1 resnorm_2 resnorm_3 resnorm
        
    % For a PID controller a second order transfer function needs to be chosen
    case {2}

        K_initial = max(y);
        tau_1_initial = t(length(t))/5;
        tau_2_initial = t(length(t))/5;
        beta_initial = t(length(t))/5;

        x_initial = [theta_initial K_initial tau_1_initial tau_2_initial];

        x = x_initial;
        
        disp('Computing model structure and parameters for PID controller.');

        [x_1, resnorm_1] = lsqnonlin(@second_order_fit_function,x,[0 min(y) 0 0],[t(length(t)) max(y) t(length(t)) t(length(t))],[],y,t);
        
        % Find dominant time constant and make it tau_1
        if x_1(3) < x_1(4)
            tau_help = x_1(3);
            x_1(3) = x_1(4);
            x_1(4) = tau_help;
            clear tau_help;
        else
        end
        
        if y(length(y)) >= y(length(y)-1)
            [x_2, resnorm_2] = lsqnonlin(@second_order_fit_function_2,[theta_initial x_1(2)/x_1(3) x_1(4)],[0 0 0],[t(length(t)) max(y) t(length(t))],[],y,t);
        else
            [x_2, resnorm_2] = lsqnonlin(@second_order_fit_function_2,[theta_initial x_1(2)/x_1(3) x_1(4)],[0 min(y) 0],[t(length(t)) 0 t(length(t))],[],y,t);
        end
        if y(length(y)) >= y(length(y)-1)
            [x_3, resnorm_3] = lsqnonlin(@second_order_fit_function_3,[theta_initial x_2(2)/x_2(3)],[0 0],[t(length(t)) max(y)],[],y,t);
        else
            [x_3, resnorm_3] = lsqnonlin(@second_order_fit_function_3,[theta_initial x_2(2)/x_2(3)],[0 min(y)],[t(length(t)) 0],[],y,t);
        end
        [x_4, resnorm_4] = lsqnonlin(@second_order_fit_function_4,[theta_initial x_1(2) x_1(3) 1.0 beta_initial],[theta_initial min(y) 0 0 0],[t(length(t)) max(y) t(length(t)) 5.0 t(length(t))],[],y,t);
        resnorm = min([resnorm_1, resnorm_2, resnorm_3, resnorm_4]);
        
        if resnorm_4 == resnorm
            % Second order transfer function plus time delay plus inverse response
            % Watch out for this option if theta * 8 < tau then process acts like an integrating process but this version does not realize it correctly
            if x_4(4) < 1.0
                K       = x_4(2);
                tau     = x_4(3);
                tau_1   = tau;
                tau_2   = 0;
                chsi    = x_4(4);
                theta   = x_4(1);
                beta    = x_4(5);
                sys_fit = tf(K*[-beta 1],[tau*tau 2*tau*chsi 1]);
                sys_fit.inputdelay = theta;
                if theta * 8 < tau
                    tau_c = max([theta, tau/(8*fac), tau_c_min, beta+theta, t(length(t))/(40.0*fac)]);
                else
                    tau_c = max([theta, tau_c_min, beta, t(length(t))/(40.0*fac)]);
                end
                func    = 9;
                integrating = 0;
            else
                resnorm = resnorm_1;
            end
        else
        end
        if resnorm_2 == resnorm
            % First order transfer function plus integrator and time delay
            K       = x_2(2);
            tau_1   = 0;
            tau_2   = x_2(3);
            theta   = x_2(1);
            beta    = 0;
            sys_fit = tf([K],conv([1 0],[tau_2 1]));
            sys_fit.inputdelay = theta;
            if theta * 8 < tau_2
                tau_c = max([theta, tau_2/(8*fac), tau_c_min, t(length(t))/(10.0*fac)]);
            else
                tau_c = max([theta, tau_c_min, t(length(t))/(10.0*fac)]);
            end
            func    = 5;
            integrating = 1;
        else
        end
        if resnorm_3 == resnorm
            % Double Integrator plus time delay
            K       = x_3(2);
            tau_1   = 0;
            tau_2   = 0;
            theta   = x_3(1);
            beta    = 0;
            sys_fit = tf([K],[1 0 0]);
            sys_fit.inputdelay = theta;
            tau_c = max([theta, tau_c_min, t(length(t))/(10.0*fac)]);
            func    = 6; 
            integrating = 1;
        else
        end
        if resnorm_1 == resnorm
            % Second order transfer function plus time delay
            K       = x_1(2);
            tau_1   = x_1(3);
            tau_2   = x_1(4);
            theta   = x_1(1);
            beta    = 0;
            sys_fit = tf([K],conv([tau_1 1],[tau_2 1]));
            sys_fit.inputdelay = theta;
            if theta * 8 < tau_1
                tau_c = max([theta, tau_1/(8*fac), tau_c_min, t(length(t))/(40.0*fac)]);
            else
                tau_c = max([theta, tau_c_min, t(length(t))/(40.0*fac)]);
            end
            func    = 2;
            integrating = 0;
        else
        end
                  
        [y_fit,t_fit] = step(sys_fit,t);
        axes(handles.axes1);
        plot(t,y,'k-',t_fit,y_fit,'k:');
        legend('Data','Second-order fit');
        xlabel('t');
        ylabel('y');
        title('y vs. t for original model and fitted model that is used for controller design');
        clear x_1 x_2 x_3 resnorm_1 resnorm_2 resnorm_3 resnorm
        
    otherwise 
        disp('Bug');
end

if data_imported == 1
    sys = sys_fit;
else
end


% Post-processing of controller structure

switch func
    case {1}

        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
        if (tau_1 > theta * 8) & ((theta * 40.0 > t(length(t))) | (theta > 0.05))
            K       = K/tau_1;
            tau_1   = 0;
            func    = 4;
        else
        end
        
    case {2}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
        if (tau_2 > theta * 8) & ((theta * 40.0 > t(length(t))) | (theta > 0.05))
            K       = K/(tau_1*tau_2);
            tau_1   = 0;
            tau_2   = 0;
            func    = 6; 
        else
        end
        
        if (tau_1 > theta * 8) & (tau_2 < theta * 8) & ((theta * 40.0 > t(length(t))) | (theta > 0.05))
            K       = K/tau_1;
            tau_1   = 0;
            func    = 5; 
        else
        end
        
    case {3}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
    case {4}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
    case {5}    
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
        if (tau_2 > theta * 8) & ((theta * 40.0 > t(length(t))) | (theta > 0.05))
            K       = K/tau_2;
            tau_2   = 0;
            func    = 6; 
        else
        end
        
    case {6}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
    case {7}
      
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
        if beta * 10.0 < tau_1
            beta    = 0;
            func    = 1;
        else
        end
        
     case {8}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end 
         
        if beta * 10.0 < max(1, theta)
            beta    = 0;
            func    = 4;
        else
        end
        
     case {9}
        
        if ((theta * 40.0 * fac < t(length(t))) | (theta < 0.01))
            theta = 0;
            tau_c = max([tau_c, t(length(t))/40.0]);
        else
        end
        
        if beta * 10.0 < tau
            beta    = 0;
        else
        end
         
     otherwise
         
       disp('Bug');  
end
 

% Print the closed-loop time constant

disp('-----------------------------------------------------------');
fprintf('Desired closed-loop time constant:  tau_c = %f\n',tau_c)         
   

% Transfer functions for each case from fitted data after post-processing
switch func
case {1}
    sys_post_fit = tf([K],[tau_1 1]);
    sys_post_fit.inputdelay = theta;
case {2}
    sys_post_fit = tf([K],conv([tau_1 1],[tau_2 1]));
    sys_post_fit.inputdelay = theta;
case {3}
    sys_post_fit = tf([K],[1]);
    sys_post_fit.inputdelay = theta;
case {4}
    sys_post_fit = tf([K],[1 0]);
    sys_post_fit.inputdelay = theta;
case {5}
    sys_post_fit = tf([K],conv([1 0],[tau_2 1]));
    sys_post_fit.inputdelay = theta;
case {6}
    sys_post_fit = tf([K],[1 0 0]);
    sys_post_fit.inputdelay = theta;
case {7}
    sys_post_fit = tf(K*[-beta 1],[tau_1 1]);
    sys_post_fit.inputdelay = theta;
case {8}
    sys_post_fit = tf(K*[-beta 1],[1 0]);
    sys_post_fit.inputdelay = theta;
case {9}
    sys_post_fit = tf(K*[-beta 1],[tau*tau 2*tau*chsi 1]);
    sys_post_fit.inputdelay = theta;
otherwise
    disp('Bug');
end;


% Compute controller and eventually retune

correctly_tuned = 0;

while correctly_tuned == 0,

    % Setting for series form PID controller
    switch func
    case {1}
        Kc      = 1/K*tau_1/(tau_c+theta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = 0;
    case {2}
        Kc      = 1/K*tau_1/(tau_c+theta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = tau_2;
    case {3}
        Kc      = 0;
        tau_i   = 0;
        tau_d   = 0;
        % This controller is a pure integral controller with no proportional action!
        Kc_tau_i= 1/(K*(tau_c+theta));
    case {4}
        Kc      = 1/(K*(tau_c+theta));
        tau_i   = 4*(tau_c+theta);
        tau_d   = 0;
    case {5}
        Kc      = 1/(K*(tau_c+theta));
        tau_i   = 4*(tau_c+theta);
        tau_d   = tau_2;
    case {6}
        Kc      = 1/(K*4*(tau_c+theta)^2);
        tau_i   = 4*(tau_c+theta);
        tau_d   = 4*(tau_c+theta);
    case {7}
        Kc      = 1/K*tau_1/(tau_c+theta+beta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = 0;
    case {8}
        Kc      = 1/(K*(tau_c+theta+beta));
        tau_i   = 4*(tau_c+theta+0.5*beta);
        tau_d   = 0;
    case {9}
        Kc      = 2*tau*chsi/(K*(tau_c+theta+beta));
        tau_i   = 2*tau*chsi;
        tau_d   = tau/(2*chsi);
    otherwise
        disp('Bug');
    end;

    % Now for ideal PI/PID controller
    if func == 3
        prop    = 0;
        integ   = Kc_tau_i;
        deriv   = 0;
        G_c     = tf([integ],[1 0]);
        G_c_non_ideal = G_c;
    else
        alpha       = 0.05;
        Kc_ideal    = Kc*(1+tau_d/tau_i);
        tau_i_ideal = tau_i*(1+tau_d/tau_i);
        tau_d_ideal = tau_d/(1+tau_d/tau_i);
        prop    = Kc_ideal;
        integ   = Kc_ideal/tau_i_ideal;
        deriv   = Kc_ideal*tau_d_ideal;
        G_c     = prop+integ*tf([1],[1 0])+deriv*tf([1 0],[1]);
        G_c_non_ideal = prop+integ*tf([1],[1 0])+Kc_ideal*tf([tau_d_ideal 1],[tau_d_ideal*alpha 1]);
    end;

    % Compute gain and phase margins to see if controller has a minimum robustness
    open_loop_sys = G_c * sys;
    warning off;
    [Gm,Pm,Wcg,Wcp] = margin(open_loop_sys);
    warning on;
    
    % If controller does not have minimum robustness then retune:
    % - increase tau_c for non-integrating systems
    % - decrease tau_c for integrating systems
    if (Gm > 1.7) & (Pm > 30.0)
        correctly_tuned = 1;
    else
        fprintf('Retuning controller:   ')
        if integrating == 0
            tau_c = tau_c * 1.2;
            fprintf('tau_c = %f\n',tau_c)
        else
            tau_c = tau_c * 0.8;
            fprintf('tau_c = %f\n',tau_c)
        end
        if (tau_c < 0.05) | (tau_c > 200.0) 
            disp('Cannot properly tune controller!');
            break;
        else
        end;
    end;        
    
end;


% Print system

disp('-----------------------------------------------------------');
disp('Original process transfer function:');
sys
disp('-----------------------------------------------------------');
disp('Fitted process transfer function:');
sys_fit
disp('-----------------------------------------------------------');
disp('Transfer function used for controller design:');
sys_post_fit
disp('-----------------------------------------------------------');
disp('Controller:');

if model_order == 1
    fprintf('G_c = P + I / s + D * s \n[P, I] = [%f, %f]\n',prop, integ)
    if func == 3
    else
        fprintf('G_c = %f*(1 + 1 / (%f * s))\n',Kc_ideal, tau_i_ideal)
    end
else    
    fprintf('G_c = P + I / s + D * s \n[P, I, D] = [%f, %f, %f]\n',prop, integ, deriv)
    fprintf('G_c = %f*(1 + 1 / (%f * s) + %f * s)\n',Kc_ideal, tau_i_ideal, tau_d_ideal)
end
disp('-----------------------------------------------------------');

set(handles.enter_tau_c,'String',num2str(tau_c));


% --------------------------------------------------------------------
function varargout = IRSetpoint_Callback(h, eventdata, handles, varargin)

global sys G_c G_c_non_ideal

closed_loop_sys = G_c * sys / (1 + G_c * pade(sys,10));
u_y_sys         = G_c_non_ideal / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = impulse(closed_loop_sys);
[u,t] = impulse(u_y_sys,t);
plot(t,y);
title('y vs. t for impulse response for the setpoint (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if (t(i) > delay) & (i ~= 1) 
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs(y(i))*(t(i+1)-t(i));
    ISE  = ISE + abs(y(i))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs(y(i))*t(i)*(t(i+1)-t(i));
end

disp('Setpoint impulse response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = SRSetpoint_Callback(h, eventdata, handles, varargin)

global sys G_c G_c_non_ideal

closed_loop_sys = G_c * sys / (1 + G_c * pade(sys,10));
u_y_sys         = G_c_non_ideal / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = step(closed_loop_sys);
[u,t] = step(u_y_sys,t);
plot(t,y);
title('y vs. t for step response for the setpoint (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if t(i) > delay
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs((1.0-y(i)))*(t(i+1)-t(i));
    ISE  = ISE + abs((1.0-y(i)))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs((1.0-y(i)))*t(i)*(t(i+1)-t(i));
end

disp('Setpoint step response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = IRLoad_Callback(h, eventdata, handles, varargin)

global sys G_c

closed_loop_sys = sys / (1 + G_c * pade(sys,10));
u_y_sys         = G_c * sys / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = impulse(closed_loop_sys);
[u,t] = impulse(u_y_sys,t);
plot(t,y);
title('y vs. t for impulse response for the load disturbance (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if (t(i) > delay) & (i ~= 1) 
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs(y(i))*(t(i+1)-t(i));
    ISE  = ISE + abs(y(i))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs(y(i))*t(i)*(t(i+1)-t(i));
end

disp('Load disturbance impulse response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = SRLoad_Callback(h, eventdata, handles, varargin)

global sys G_c

closed_loop_sys = sys / (1 + G_c * pade(sys,10));
u_y_sys         = G_c * sys / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = step(closed_loop_sys);
[u,t] = step(u_y_sys,t);
plot(t,y);
title('y vs. t for step response for the load disturbance (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if t(i) > delay
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs(y(i))*(t(i+1)-t(i));
    ISE  = ISE + abs(y(i))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs(y(i))*t(i)*(t(i+1)-t(i));
end

disp('Load disturbance step response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = IRDist_Callback(h, eventdata, handles, varargin)

global sys G_c G_c_non_ideal

closed_loop_sys = 1 / (1 + G_c * pade(sys,10));
u_y_sys         = G_c_non_ideal / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = impulse(closed_loop_sys);
[u,t] = impulse(u_y_sys,t);
plot(t,y);
title('y vs. t for impulse response for the output disturbance (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if (t(i) > delay) & (i ~= 1) 
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs(y(i))*(t(i+1)-t(i));
    ISE  = ISE + abs(y(i))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs(y(i))*t(i)*(t(i+1)-t(i));
end

disp('Output disturbance impulse response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = SRDist_Callback(h, eventdata, handles, varargin)

global sys G_c G_c_non_ideal

closed_loop_sys = 1 / (1 + G_c * pade(sys,10));
u_y_sys         = G_c_non_ideal / (1 + G_c * pade(sys,10));
axes(handles.axes1);
[y,t] = step(closed_loop_sys);
[u,t] = step(u_y_sys,t);
plot(t,y);
title('y vs. t for step response for the output disturbance (Pade approximation is used)');

IAE  = 0;
ISE  = 0;
ITAE = 0;
TV   = 0;
delay= get(sys,'InputDelay');
for i = 1:length(y)-1
    if t(i) > delay
        TV   = TV + abs(u(i+1) - u(i));
    else 
    end
    IAE  = IAE + abs(y(i))*(t(i+1)-t(i));
    ISE  = ISE + abs(y(i))^2*(t(i+1)-t(i));
    ITAE = ITAE + abs(y(i))*t(i)*(t(i+1)-t(i));
end

disp('Output disturbance step response:');
fprintf('IAE  = %f\nISE  = %f\nITAE = %f\nTV   = %f\n',IAE, ISE, ITAE, TV)
disp('-----------------------------------------------------------');

% --------------------------------------------------------------------
function varargout = Close_Callback(h, eventdata, handles, varargin)

close;


% --------------------------------------------------------------------
function varargout = Update_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys Kc tau_c tau_i tau_d prop integ deriv K theta tau_1 tau_2 integrating G_c_non_ideal tau chsi beta

% Compute controller and eventually retune

correctly_tuned = 0;

while correctly_tuned == 0,

    % Setting for series form PID controller
    switch func
    case {1}
        Kc      = 1/K*tau_1/(tau_c+theta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = 0;
    case {2}
        Kc      = 1/K*tau_1/(tau_c+theta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = tau_2;
    case {3}
        Kc      = 0;
        tau_i   = 0;
        tau_d   = 0;
        % This controller is a pure integral controller with no proportional action!
        Kc_tau_i= 1/(K*(tau_c+theta));
    case {4}
        Kc      = 1/(K*(tau_c+theta));
        tau_i   = 4*(tau_c+theta);
        tau_d   = 0;
    case {5}
        Kc      = 1/(K*(tau_c+theta));
        tau_i   = 4*(tau_c+theta);
        tau_d   = tau_2;
    case {6}
        Kc      = 1/(K*4*(tau_c+theta)^2);
        tau_i   = 4*(tau_c+theta);
        tau_d   = 4*(tau_c+theta);
    case {7}
        Kc      = 1/K*tau_1/(tau_c+theta+beta);
        tau_i   = min(tau_1, 4*(tau_c+theta));
        tau_d   = 0;
    case {8}
        Kc      = 1/(K*(tau_c+theta+beta));
        tau_i   = 4*(tau_c+theta+0.5*beta);
        tau_d   = 0;
    case {9}
        Kc      = 2*tau*chsi/(K*(tau_c+theta+beta));
        tau_i   = 2*tau*chsi;
        tau_d   = tau/(2*chsi);
    otherwise
        disp('Bug');
    end;

    % Now for ideal PI/PID controller
    if func == 3
        prop    = 0;
        integ   = Kc_tau_i;
        deriv   = 0;
        G_c     = tf([integ],[1 0]);
        G_c_non_ideal = G_c;
    else
        alpha       = 0.05;
        Kc_ideal    = Kc*(1+tau_d/tau_i);
        tau_i_ideal = tau_i*(1+tau_d/tau_i);
        tau_d_ideal = tau_d/(1+tau_d/tau_i);
        prop    = Kc_ideal;
        integ   = Kc_ideal/tau_i_ideal;
        deriv   = Kc_ideal*tau_d_ideal;
        G_c     = prop+integ*tf([1],[1 0])+deriv*tf([1 0],[1]);
        G_c_non_ideal = prop+integ*tf([1],[1 0])+Kc_ideal*tf([tau_d_ideal 1],[tau_d_ideal*alpha 1]);
    end;

    % Compute gain and phase margins to see if controller has a minimum robustness
    open_loop_sys = G_c * sys;
    warning off;
    [Gm,Pm,Wcg,Wcp] = margin(open_loop_sys);
    warning on;
    
    % If controller does not have minimum robustness then retune:
    % - increase tau_c for non-integrating systems
    % - decrease tau_c for integrating systems
    if (Gm > 1.05) & (Pm > 5.0)
        correctly_tuned = 1;
    else
        fprintf('Retuning controller:   ')
        if integrating == 0
            tau_c = tau_c * 1.2;
            fprintf('tau_c = %f\n',tau_c)
        else
            tau_c = tau_c * 0.8;
            fprintf('tau_c = %f\n',tau_c)
        end
        if (tau_c < 0.05) | (tau_c > 200.0) 
            disp('Cannot properly tune controller!');
            break;
        else
        end;
    end;        
    
end;

fprintf('Retuning: \nDesired closed-loop time constant:  tau_c = %f\n',tau_c)  

disp('Controller:');

if model_order == 1
    fprintf('G_c = P + I / s + D * s \n[P, I] = [%f, %f]\n',prop, integ)
    if func == 3
    else
        fprintf('G_c = %f*(1 + 1 / (%f * s))\n',Kc_ideal, tau_i_ideal)
    end
else    
    fprintf('G_c = P + I / s + D * s \n[P, I, D] = [%f, %f, %f]\n',prop, integ, deriv)
    fprintf('G_c = %f*(1 + 1 / (%f * s) + %f * s)\n',Kc_ideal, tau_i_ideal, tau_d_ideal)
end
disp('-----------------------------------------------------------');

set(handles.enter_tau_c,'String',num2str(tau_c));

% --------------------------------------------------------------------
function varargout = enter_tau_c_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys Kc tau_c tau_i tau_d prop integ deriv

% Get the new value for tau_c
NewStrVal = get(handles.enter_tau_c,'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0.05) | (NewVal>100),
    % Revert to last value, as indicated by Slider
    set(handles.enter_tau_c,'String',num2str(tau_c));
else
    % Update tau_c
    tau_c = NewVal;
end