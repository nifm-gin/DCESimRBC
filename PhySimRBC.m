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

function [C, Cart, Cmax, Tstart, T_base] = PhySimRBC(Model,Seq,t,too,G,C,Scale_factor)

% Compute the evolution of the physiology/physic related lattices
% 
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012


%% Use a linear increasing CA concentration
T_base = 5;
Cmax = 18/Scale_factor;
Tstart = T_base+Seq.Tacq;
if any(t <= Tstart), Cart(t <= Tstart) = 0; end
if any(t > Tstart), Cart(t > Tstart) = (Cmax/(Model.Tmax-Tstart))*(t(t > Tstart)-Tstart); end
C(logical(ones(size(G.vasc.P)) - full(G.vasc.P))) = Cart;
Cv = 0;

%% Gd Permeability in EES
C = C + (Model.phy.CA.k * too.PGwghtI2O .* (Cv - C)*Model.dt) .* full(G.vasc.Pp);
    
%% Gd Diffusion
C = (abs(ifftn(fftn(C.*full(G.ees.P)) .* too.FTDGker)) + C.*full(G.ees.P).*too.DGwgh) .* full(G.ees.P) + full(G.vasc.P)*Cv;
