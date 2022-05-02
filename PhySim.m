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

function C = PhySim(Model,t,too,G,C)

% Compute the evolution of the physiology/physic related lattices
% 
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012


%% Compute AIF
ca = AIFGen(Model.AIF,t,Model.dt);

if Model.phy.CA.k ~= 0 % vessels are permeable
    %% Vacular
    if too.Fdt < too.vf
        Cv    = mean( C( logical( full( G.vasc.P ) ) ) );
        Cperi = C(logical(full(G.vasc.Pp)));
        dCperi = Model.phy.CA.k*too.PGwghtI2O(logical(full(G.vasc.Pp))) .* ( Cv - Cperi);
        dCp = sum(dCperi(:))/sum(G.vasc.P(:));
        Cv = Cv + too.Fdt/too.vf * (ca - Cv) - dCp*Model.dt;
    else
        Cv = ca;
    end
    C(logical(full(G.vasc.P))) = Cv;
    
    %% Gd Permeability in EES
    C = C + (Model.phy.CA.k * too.PGwghtI2O .* (Cv - C)*Model.dt) .* full(G.vasc.Pp);
    
    %% Gd Diffusion
    C = (abs(ifftn(fftn(C.*full(G.ees.P)) .* too.FTDGker)) + C.*full(G.ees.P).*too.DGwgh) .* full(G.ees.P) + full(G.vasc.P)*Cv;
else
    %% Vacular
    if ca~=0
        '';
    end
    if too.Fdt < too.vf
        Cv = mean(C(logical(full(G.vasc.P))));
        Cv = Cv + too.Fdt/too.vf * (ca - Cv);
    else
        Cv = ca;
    end
    C(logical(full(G.vasc.P))) = Cv;
end
