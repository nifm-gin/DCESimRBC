% Copyright 2012, 2013 Nicolas Pannetier, Clément Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
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

function [Model, Seq] = ReadModelAndSeq(File)

% Read and parse the parameter file
% Input: 
%   - File  : Path to the parameter file
% Output:
%   - Model     : Structure containing the voxel related parameters
%   - Seq       : Structure containing the sequence related parameters
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

Model = struct;
Seq = struct;

%% read lines
fid = fopen(File,'rt');
C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
fclose(fid);

%% Search for Model and Seq structure in the file
for a=1:numel(C)
    ModelLine(a) = strncmp(C{a},'Model.',6);
    SeqLine(a) = strncmp(C{a},'Seq.',4);
end
MLine = find(ModelLine);
SLine = find(SeqLine);

%% Read Model strucutre
for a=1:numel(MLine)
   eval(C{MLine(a)}); 
end

%% Read Seq structure
for a=1:numel(SLine)
   eval(C{SLine(a)}); 
end