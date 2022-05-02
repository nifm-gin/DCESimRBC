% Copyright 2013 Nicolas Pannetier, Clément Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
% This file is part of DCESIM.
% 
%   DCESIM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%    DCESIM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DCESIM.  If not, see <http://www.gnu.org/licenses/>.

function [progress_str, tic_dt] = DisplayProgress(tt,Model,progress_str,tic_dt)

% Display algorithm progression
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

percentFinished = floor(tt/Model.Tmax*100);
elapsed_time = toc;

if isempty(progress_str)
    prevLength = 0;
    time_left = 0;
    rate = 0;
else
    prevLength = numel(progress_str);
    time_left     = (Model.Tmax/tt-1) * elapsed_time;
    elapsed_time_dt = toc(tic_dt);
    rate = elapsed_time_dt;
end
progress_str = sprintf('%3.0f %% simulated in %3.0f s, at %3.2f s/dt (%3.0f /%3.0f ms done)\n Estimated time left: %3.0f s \n',...
    percentFinished,elapsed_time,rate,tt*1000,Model.Tmax*1000,time_left);
fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
tic_dt=tic;