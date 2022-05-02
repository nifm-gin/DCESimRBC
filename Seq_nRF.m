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

function [M, Seqflag,too] = Seq_nRF(tt,dt,para,M,too)

% Simulate an n RF pulse sequence
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% NRF with read cst gradient
Tacq = para.Tacq;
TR = para.TR;
TRF = para.RF.exc.time; 

%% RF pulse
if any( abs(TRF - mod(tt,TR))< dt *(1-1e-9))
    
    ind = find(abs(TRF - mod(tt,TR))< dt *(1-1e-9));
        
    Mtmp = reshape(cat(3,real(M.per(:)),imag(M.per(:)),M.par(:)),[numel(M.per(:)) 3]);
    Mtmp = [too.RFtrans(1,ind) .* Mtmp(:,1) + too.RFtrans(2,ind) .* Mtmp(:,2) + too.RFtrans(3,ind) .* Mtmp(:,3), ...
        too.RFtrans(4,ind) .* Mtmp(:,1) + too.RFtrans(5,ind) .* Mtmp(:,2) + too.RFtrans(6,ind) .* Mtmp(:,3), ...
        too.RFtrans(7,ind) .* Mtmp(:,1) + too.RFtrans(8,ind) .* Mtmp(:,2) + too.RFtrans(9,ind) .* Mtmp(:,3)];
    Mtmp = reshape(Mtmp,[size(M.per,1) size(M.per,2) 3]);
    
    M.par = Mtmp(:,:,3);
    M.per = Mtmp(:,:,1) + 1i * Mtmp(:,:,2);
    too.inc = too.inc+1;
end

%% Test pour acquisition
if any( abs(Tacq - mod(tt,TR))< dt *(1-1e-9))
    Seqflag.Acq = 1;
else
    Seqflag.Acq = 0;
end