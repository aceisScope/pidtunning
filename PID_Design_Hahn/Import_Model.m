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


function varargout = Import_Model(varargin)
% IMPORT_MODEL Application M-file for Import_Model.fig
%    FIG = IMPORT_MODEL launch Import_Model GUI.
%    IMPORT_MODEL('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 23-Sep-2002 20:52:42

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    % Populate the listbox
    %set(handles.y,'Value',[])
    %set(handles.t,'Value',[])
    update_listbox(handles)
    
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
function varargout = num_Callback(h, eventdata, handles, varargin)

global y t sys


% --------------------------------------------------------------------
function varargout = den_Callback(h, eventdata, handles, varargin)

global y t sys


% --------------------------------------------------------------------
function varargout = Import_Callback(h, eventdata, handles, varargin)

global y t sys data_imported

y = 0;
t = 0;

num = eval(get(handles.num,'String'));
den = eval(get(handles.den,'String'));
delay = eval(get(handles.delay,'String'));

sys = tf(num,den);
sys.inputdelay = delay;
[y,t] = step(sys);

% Make sure that t = delay is one of the outputs of the data
% This is still experimental since it seems to cause numerical problems lateron for the optimization
% Conclusion from this is that somehow the time delay is causing serious problems for the integrator
% Maybe one should rewrite the problems such that time delay is set as best as possible and the rest is simply fitted?

%i = 1;
%while t(i) < delay,
%    if i < length(t)
%        i = i + 1;
%    else
%        break;
%    end
%end
%
%t = [t(1:i-1);delay;t(i:length(t))];
%y = [y(1:i-1);y(i-1);y(i:length(y))];

data_imported = 0;

disp('Model imported.');
disp('-----------------------------------------------------------');

close;


% --------------------------------------------------------------------
function varargout = delay_Callback(h, eventdata, handles, varargin)

global y t sys

% --------------------------------------------------------------------
function varargout = t_Callback(h, eventdata, handles, varargin)

global y t sys


% --------------------------------------------------------------------
function varargout = y_Callback(h, eventdata, handles, varargin)

global y t sys


% --------------------------------------------------------------------
function varargout = Import2_Callback(h, eventdata, handles, varargin)

global y t sys data_imported error

y = 0;
t = 0;
error = 0;

[y,t] = get_var_names(handles);

if error == 1
    data_imported = 0;
else
    data_imported = 1;

    disp('Data for model imported.');
    disp('-----------------------------------------------------------');

    clear error
    close
end


function update_listbox(handles)
% Updates the listbox to match the current workspace
  vars = evalin('base','who');
  set(handles.t,'String',vars)
  set(handles.y,'String',vars)
  
function [y,t] = get_var_names(handles)

global error

% Returns the names of the two variables to plot
list_entries_t = get(handles.t,'String');
index_selected_t = get(handles.t,'Value');
list_entries_y = get(handles.y,'String');
index_selected_y = get(handles.y,'Value');
if length(index_selected_t) ~= 1
	errordlg('You must select one variable','Incorrect Selection','modal')
    error = 1.0;
else
	t = evalin('base',list_entries_t{index_selected_t});
end 
[n_t,m_t] = size(t);
if (m_t ~= 1) & (n_t ~= 1)
    disp('Wrong dimension for imported variable t.');
    disp('-----------------------------------------------------------');
    error = 1.0;
else
end
if (m_t ~= 1) & (n_t == 1)
    t = t';
else
end
if ~isequal(t,sort(t)) & (error ~= 1.0)
    disp('t must be in ascending order.');
    disp('-----------------------------------------------------------');
    error = 1.0;
else
end
if (min(t) < 0) & (error ~= 1.0) 
    disp('t cannot be negative.');
    disp('-----------------------------------------------------------');
    error = 1.0;
else
end
if length(index_selected_y) ~= 1
	errordlg('You must select one variable','Incorrect Selection','modal')
    error = 1.0;
else
	y = evalin('base',list_entries_y{index_selected_y});
end 
[n_y,m_y] = size(y);
if (m_y ~= 1) & (n_y ~= 1)
    disp('Wrong dimension for imported variable y.');
    disp('-----------------------------------------------------------');
    error = 1.0;
else
end
if (m_y ~= 1) & (n_y == 1)
    y = y';
else
end
if (length(y) ~= length(t)) & (error ~= 1.0)
    disp('y and t must have same length.');
    disp('-----------------------------------------------------------');
    error = 1.0;
else
end