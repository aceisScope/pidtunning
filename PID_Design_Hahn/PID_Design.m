% PI/PID controller tuning 
% Copyright (C) 2002: Juergen Hahn
% Department of Chemical Engineering
% Texas A&M University
% College Station, TX 77843-3122
% USA
% hahn@che.tamu.edu
%
% This file is intended for use with MATLAB 6.1 or later versions. The
% optimization toolbox is also required.
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


function varargout = PID_Design(varargin)
% PID_DESIGN Application M-file for PID_Design.fig
%    FIG = PID_DESIGN launch PID_Design GUI.
%    PID_DESIGN('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 16-Sep-2002 14:49:10

global y t model_order func prop integ deriv G_c sys_fit sys Kc tau_c tau_i tau_d prop integ deriv K theta tau_1 tau_2 beta integrating G_c_non_ideal data_imported

if nargin == 0
    clc;  
    disp('-----------------------------------------------------------');
    clear global
else
end   
        
        
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

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
function varargout = Design_Controller_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys

Controller_Design;



% --------------------------------------------------------------------
function varargout = Analyze_Controller_Callback(h, eventdata, handles, varargin)

global y t model_order func prop integ deriv G_c sys_fit sys

analyze_controller;


% --------------------------------------------------------------------
function varargout = Import_Model_Callback(h, eventdata, handles, varargin)

global y t sys

Import_Model;


% --------------------------------------------------------------------
function varargout = Quit_Callback(h, eventdata, handles, varargin)

close;



% --------------------------------------------------------------------
function varargout = About_Callback(h, eventdata, handles, varargin)

about;