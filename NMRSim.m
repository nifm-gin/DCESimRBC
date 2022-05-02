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

function [M, Seqflag, save] = NMRSim(Model,too,G,C,M,Pulseq,Seq,tt)

% Compute the evolution of the magnetization matrices
% 
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012



%% Compute relevant matrices

% khi
khi =   Model.phy.vasc.khi * full(G.vasc.P) + Model.phy.CA.khi_mol * C .* full(G.vasc.P) + ...
        Model.phy.cell.khi * full(G.cell.P) + Model.phy.CA.khi_mol * C .* full(G.cell.P) + ...
        Model.phy.ees.khi  * full(G.ees.P)  + Model.phy.CA.khi_mol * C .* full(G.ees.P);
    
% R2
R2 =    1/Model.phy.vasc.T2 * full(G.vasc.P) + Model.phy.CA.r2 * C.* full(G.vasc.P) + ...
        1/Model.phy.cell.T2 * full(G.cell.P) + Model.phy.CA.r2 * C.* full(G.cell.P) + ...
        1/Model.phy.ees.T2  * full(G.ees.P)  + Model.phy.CA.r2 * C.* full(G.ees.P);
% R1    
R1 =    1/Model.phy.vasc.T1 * full(G.vasc.P) + Model.phy.CA.r1 * C.* full(G.vasc.P) + ...
        1/Model.phy.cell.T1 * full(G.cell.P) + Model.phy.CA.r1 * C.* full(G.cell.P) + ...
        1/Model.phy.ees.T1  * full(G.ees.P)  + Model.phy.CA.r1 * C.* full(G.ees.P);

% Compute magnetic field perturbations dB
if Model.Flag.B0Ori3D == 1
    thetaB0 = [pi/2 pi/2 0];
    [b1, m1, kk1] = unique(thetaB0, 'first');
    for p=1:numel(b1)
        warning('off');
        FtB(:,:) = Model.phy.B0 * fftshift(fftn(khi)) .*(1/3 - ((too.yB * sin(thetaB0(m1(p)))).^2) ./ (too.xB.^2+too.yB.^2));
        warning('on');
        %     FTB(x_pix/2+1,y_pix/2+1) = 0;                      %Pour normalisation -B0*xe/3 ou 0
        B = real(ifftn(ifftshift(FtB)));
        dB(:,:,p)=B;
    end
    dBm = mean(dB(:,:,kk1),3);
else
    warning('off');
    FtB = phy.B0 * fftshift(fftn(khi)) .*(1/3 - ((too.yB * sin(pi/2)).^2) ./ (too.xB.^2+too.yB.^2));
    warning('on');
%     FtB(size(FtB,1)/2+1,size(FtB,2)/2+1) = 0;                      %Pour normalisation -B0*xe/3 ou 0
    B = real(ifftn(ifftshift(FtB)));
    dBm=B;
end

% Pool and compute rotation matrix in transverse plan
rot = exp( (-1i* Model.phy.gamma *dBm - R2) * Model.dt);

%% Evolution
M.per = M.per.*rot;
M.par = (M.par-M.M0) .* exp(-Model.dt*R1) + M.M0; 

%% Compute diffusion of water
if Model.phy.DH2O ~= 0
    M.par = ifftn(fftn(M.par) .* too.FTDHker);
    M.per = ifftn(fftn(M.per) .* too.FTDHker);
end

%% Play MR sequence 
[M, Seqflag] = Pulseq(tt,Model.dt,Seq,M,too);



