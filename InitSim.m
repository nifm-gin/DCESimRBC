% Copyright 2013 Nicolas Pannetier, Cl�ment Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
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

% Copyright 2012, 2013 Nicolas Pannetier, Cl�ment Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
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

function [too, C, M] = InitSim( Model , G)

% Initiate the matrices and structure used in the simulation
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012


%% Initialize variables
if Model.Project == 1
    fov = 70e-06;       %Fixed FOV for RBC project
else
    fov = Model.geo.vasc.R * sqrt( Model.geo.vasc.N * pi / Model.geo.vasc.vp);
end

dx = fov / Model.geo.res;

[x,y] = ndgrid( -fov/2 : dx : fov/2-dx , -fov/2:dx:fov/2-dx);
too.x = x;
too.y = y;
         
[xB,yB] = ndgrid( -fov/2 + dx/2 : dx : fov/2-dx/2,-fov/2+dx/2:dx:fov/2-dx/2);
too.xB = xB;
too.yB = yB;
   

%% Gd Diffusion kernel
sig = sqrt( 2 * Model.phy.CA.D * Model.dt );
too.DGker = 1/(2*pi*sig^2)^(2/2) * exp(- (x.^2 + y.^2)/(2*sig^2));
too.DGker = too.DGker/sum(too.DGker(:));
too.FTDGker = fftn(fftshift(too.DGker));

%% Gd diffusion weighting matrix
too.DGwgh = abs(ifftshift(ifftn(fftn(full(G.vasc.P+G.cell.P)).*fftn(too.DGker)))).*full(G.ees.P);

%% Weighting matrix at vessel/interstitium interface
% from in to out
S = conv2(full(G.vasc.P),[0 1 0; 1 1 1; 0 1 0],'same') .* (ones(size(G.vasc.P)) - full(G.vasc.P));
too.PGwghtI2O = sum(full(G.ees.P(:))) * S/sum(S(:));

% from out to in
S = conv2(full(1- G.vasc.P),[0 1 0; 1 1 1; 0 1 0],'same') .* full(G.vasc.P);
too.PGwghtO2I = sum(full(G.ees.P(:))) * S/sum(S(:));

%% Water diffusion kernel
sig = sqrt(2*Model.phy.DH2O*Model.dt);
too.DHker = 1/(2*pi*sig^2)^(2/2) * exp(- (x.^2 + y.^2)/(2*sig^2));
too.DHker = too.DHker/sum(too.DHker(:));
too.FTDHker = fftn(fftshift(too.DHker));

%% Initialization of matrices
C = zeros(size(G.vasc.P));
M.M0 = Model.phy.M0.vasc * full(G.vasc.P) + Model.phy.M0.ees * full(G.cell.P) + Model.phy.M0.ees * full(G.ees.P);
M.M0 = M.M0/sum(M.M0(:));
M.par = M.M0;
M.per = zeros(size(G.vasc.P));

too.vf = Model.geo.vasc.vp;
too.Fdt = Model.AIF.F*Model.dt;
too.gamma = Model.phy.gamma;
too.dt = Model.dt;
too.inc = 1;

