% Copyright 2012, 2013 Nicolas Pannetier, Clement Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
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

function [M, Seqflag,too] = Seq_RFspoil(tt,dt,para,M,too)
%
% Simulate an n RF pulse sequence
%
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% NRF with read cst gradient
Tacq = para.Tacq;
TR = para.TR;
TRF = para.RF.exc.time;
ang = para.RF.exc.ang*pi/180;
pha = ((too.inc * (too.inc-1)/2) * para.Spoil.RF.phainc + para.RF.exc.pha)*pi/180;

if ~isfield(too,'Spoil')
    too.Spoil.inc = 1;
    too.Spoil.phi = 0;
end

%% RF pulse
if any( abs(TRF - mod(tt,TR))< dt *(1-1e-9))
      
    ra =[cos(pha) sin(pha) 0 ];
    rg = ang;

    RFtrans = zeros(3,3);
    RFtrans(1) = ra(1).^2 + (1 - ra(1).^2) .* cos (rg);
    RFtrans(2) = ra(1).*ra(2) .* (1 -cos(rg)) - ra(3) .* sin(rg);
    RFtrans(3) = ra(1).*ra(3) .* (1 -cos(rg)) + ra(2) .* sin(rg);
    RFtrans(4) = ra(1).*ra(2) .* (1 -cos(rg)) + ra(3) .* sin(rg);
    RFtrans(5) = ra(2).^2 + (1 - ra(2).^2) .* cos (rg);
    RFtrans(6) = ra(2).*ra(3) .* (1 -cos(rg)) - ra(1) .* sin(rg);
    RFtrans(7) = ra(1).*ra(3) .* (1 -cos(rg)) - ra(2) .* sin(rg);
    RFtrans(8) = ra(2).*ra(3) .* (1 -cos(rg)) + ra(1) .* sin(rg); 
    RFtrans(9) = ra(3).^2 + (1 - ra(3).^2) .* cos (rg);
    
    Mtmp = cat(1, real( M.per( : ) ).', imag( M.per( : ) ).', M.par( : ).' );
    Mtmp = RFtrans.' * Mtmp;
    Mtmp = reshape(Mtmp,[3, size(M.per) ]);
    
    M.par = permute( Mtmp(3,:,:,:), [ 2 3 4 1]);
    M.per = permute( complex( Mtmp(1,:,:,:) , Mtmp(2,:,:,:)), [2 3 4 1]);
    too.inc = too.inc+1;
end

%% Grad Spoiling
dspoil = para.Spoil.Grad.dur;
if (TR - mod(tt,TR))  <= dspoil*(1+1e-9)
        M.per = M.per.*too.SpoilGradrot;
end

%% Test pour acquisition
if any( abs(Tacq - mod(tt,TR))< dt *(1-1e-9))
    Seqflag.Acq = 1;
else
    Seqflag.Acq = 0;
end