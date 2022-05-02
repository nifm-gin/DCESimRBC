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

function [too] = InitSeq(para,too)

% Compute RF pulses rotation operator
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

ang = para.RF.exc.ang*pi/180;
pha = para.RF.exc.pha*pi/180;    

for ind=1:numel(para.RF.exc.ang)
    
    ra =[cos(pha(ind)) sin(pha(ind)) 0 ];
    rg = ang(ind);

    too.RFtrans(1,ind) = ra(1).^2 + (1 - ra(1).^2) .* cos (rg);
    too.RFtrans(2,ind) = ra(1).*ra(2) .* (1 -cos(rg)) - ra(3) .* sin(rg);
    too.RFtrans(3,ind) = ra(1).*ra(3) .* (1 -cos(rg)) + ra(2) .* sin(rg);
    too.RFtrans(4,ind) = ra(1).*ra(2) .* (1 -cos(rg)) + ra(3) .* sin(rg);
    too.RFtrans(5,ind) = ra(2).^2 + (1 - ra(2).^2) .* cos (rg);
    too.RFtrans(6,ind) = ra(2).*ra(3) .* (1 -cos(rg)) - ra(1) .* sin(rg);
    too.RFtrans(7,ind) = ra(1).*ra(3) .* (1 -cos(rg)) - ra(2) .* sin(rg);
    too.RFtrans(8,ind) = ra(2).*ra(3) .* (1 -cos(rg)) + ra(1) .* sin(rg); 
    too.RFtrans(9,ind) = ra(3).^2 + (1 - ra(3).^2) .* cos (rg);
end

%% Grad Spoiling
if isfield(para,'Spoil')
   maxX = max(too.xB(:));
   dxX = diff(too.xB(1:2));
   G = para.Spoil.Grad.A*2*pi/(2*maxX+dxX);
   too.SpoilGradrot = exp( - 1i* G * (too.xB + maxX));
end
