% PI/PID controller tuning 
% Copyright (C) 2002: Juergen Hahn
% Department of Chemical Engineering
% Texas A&M University
% College Station, TX 77843-3122
% USA
% hahn@che.tamu.edu
%
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
%
% This file is intended for use with MATLAB 6.1 or later versions. The
% optimization toolbox is also required.
%
%
% Controllers of PID type are by far the most widely used controllers 
% in the chemical process industries. However, in many cases the 
% controllers are not properly tuned resulting in a significant loss of 
% performance. In order to counter this, new and simple tuning rules have 
% been derived for single-loop PID controllers. These rules can be automated 
% in order to assist with the design of PI and PID controllers. This 
% software implements existing PI and PID tuning rules for processes that 
% are either given in transfer function form or where data for the 
% input-output behavior of the process is available. The controller tuning 
% and evaluation of the closed-loop performance is performed interactively 
% using a graphical user interface in order to enhance ease of use.
%
% The tuning algorithm is based upon a modification of Skogestad's SIMC 
% routines. One change is that the transfer functions are simplified 
% based upon numerically fitting a function to the data that are 
% generated from the original model. This way it is also possible to 
% tune a controller if only data is available, but no transfer function 
% representation has been determined. Additionally, tuning methods 
% for underdamped second order systems and systems which exhibit 
% inverse response have been added.
