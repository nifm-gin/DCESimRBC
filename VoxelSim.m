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

function [S, C, t, Cblood_inc, Relax2_inc,Cblood,Relax2,Scale_factor] = VoxelSim(File)

% This code simulates the MR signal of DCE-like experiments
%
% Input:
%   - File  : Path to the parameter file
% Output:
%   - S     : MR Signal
%   - C     : CA concentration
%   - t     : Acquisition time
%
% Original DCE code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% Adjusted DSC code for intra-arterial experiments by D. van Dorth, Leiden University Medical Center, 2021
% Abstract publication at ESMRM-B 2020: 'Relationship between R2(*) Relaxation and Gd-DTPA Concentration in Whole Blood Assessed by Simulations'; D. van Dorth, L. Hirschler, K. Venugopal, B. Schmitz Abecassis, D. H. J. Poot, M. Smits, J. A. Hernandez-Tamames, M. J. P. van Osch; ESMRM-B 2020; abstract number A-1229
% Abstract publication at ISMRM 2021: 'Dependency of R2(*) Relaxation on Gd-DTPA Concentration in Arterial Blood: Influence of Hematocrit and Magnetic Field Strength'; D. van Dorth, K. Venugopal, D. H. J. Poot, M. Smits, J. H. J. M. de Bresser, J. A. Hernandez-Tamames, M. J. P. van Osch; ISMRM 2021 

%%  INIT Simulation
%%###################
%% Load Input Parameters
[Model, Seq] = ReadModelAndSeq(File);

%% Define random number seed
stream = RandStream('mt19937ar');
reset(stream,Model.geo.vasc.Id*1e6+Model.geo.cell.Id);

%% NMR SEQUENCE
Pulseq = eval(Seq.Name); % create function handle. 

%% Other paramters
Model.phy.gamma = 2.675*10^8;
if Model.Project == 1       %In case of RBC project, use a fixed RBC susceptibility khi0
    Model.phy.vasc.khi = Model.phy.vasc.khi0;
else
    Model.phy.vasc.khi = 4*pi*Model.phy.vasc.khi0*Model.phy.vasc.Hct*(1-Model.phy.vasc.Y);
end
tic_dt = 0;
tic

%%  START Simulation
%%###################

%% Generate Geometry lattice
G = GeoGen(Model);

%% Generate static matrices
% For the voxel
[too, CA, M] = InitSim(Model,G);
% For the MR sequence
too = InitSeq(Seq,too);

%% Other Parameters
count = 1;
C = [];
S = [];
t = [];
str = [];
AlreadyDraw = 0;

%% Extra parameters in case of RBC project
if Model.Project == 1
    fov = 70e-6;
    Cblood = [];
    Relax2 = [];
    if Model.Shape == 0
        HCT = Model.geo.vasc.N*pi*Model.geo.vasc.R^2/fov^2; %Calculate hematocrit value for circular RBCs
    else
        HCT = Model.geo.vasc.N*pi*Model.geo.vasc.R*Model.geo.vasc.R2/fov^2; %Calculate hematocrit value for ellipsoidal RBCs
    end
    Scale_factor = 1-HCT;
end

%% Main Time Loop
for tt=0:Model.dt:Model.Tmax
    
    %% Print progression
    if Model.Flag.verbose
        [str, tic_dt] = DisplayProgress(tt,Model,str,tic_dt);
    end
    
    %% Display real time lattice M, [CA] and S in case of vascular project
    if Model.Project == 0
        if Model.Flag.display
            if ~AlreadyDraw
                FigProp = DisplayFig_gen(G,M,CA,S,C,t);
                AlreadyDraw = 1;
            else
                FigProp = DisplayFig_up(G,M,CA,S,C,t,FigProp,tt);
            end
        end
    end
    
    %% Physiology block
    if Model.Project == 1 
        if tt-fix(tt)<=Seq.Tacq
            [CA, CB, Cmax, Tstart, T_base] = PhySimRBC(Model,Seq,tt,too,G,CA,Scale_factor);
        end
    else
        CA = PhySim(Model,tt,too,G,CA);
    end
    
    %% NMR block
    [M, Seqflag] = NMRSim(Model,too,G,CA,M,Pulseq,Seq,tt);
    
    %% Acquisition
    if tt == 0
        max_mag = max(abs(M.per(:)));
        max_phase = abs(max(angle(M.per(:))));
    end
    
    if Seqflag.Acq
        t(count) = tt;
        S(count) = sum( M.per( : ) );
        C(count) = mean( CA( : ) );
        if Model.Project == 1
            Cblood(count) = CB;
            if tt>Tstart
                S_base = mean(S(3:4));      %Calculate the baseline signal
                Relax2(count) = -(1/Seq.Tacq)*log(abs(S(count))/abs(S_base));   %Calculate the R2* relaxation from the MR signal magnitude
            end
        end
        count = count+1;
        
        %% RBC project: Plotting MR signal distribution and CA concentration at specified time  
        if  Model.Project == 1 && tt == Model.Tmax
            figure,
            VessContour = contourc(full(double(G.vasc.P)),1);
            plot(0,0),
            h1 = imagesc(abs(M.per)); hold on; axis square, axis off,title('Magnitude'), caxis([0,max_mag]), colorbar %;colormap gray
            inds = 1;
            hold on
            while inds < size(VessContour,2)
                plot(0,0),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
                inds = inds + VessContour(2,inds)+1;
            end
            hold off
            
            figure,
            plot(0,0),
            h2 = imagesc(angle(M.per)); hold on; axis square, axis off,title('Phase'), caxis([-max_phase,max_phase]), colorbar
            inds = 1;
            hold on,
            while inds < size(VessContour,2)
                plot(0,0),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
                inds = inds + VessContour(2,inds)+1;
            end
            hold off
            
            figure,
            plot(0,0),
            tmp = CA;
            h3 = imagesc(tmp); hold on; axis square, axis off,title('[CA] (mM)'), caxis([0 Cmax]), colorbar
            inds = 1;
            hold on,
            while inds < size(VessContour,2)
                plot(0,0),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
                inds = inds + VessContour(2,inds)+1;
            end
            hold off
            
            FigProp.h1 =h1;
            FigProp.h2 =h2;
            FigProp.h3 =h3;
            
            set(FigProp.h1,'CData',abs(M.per))
            set(FigProp.h2,'CData',angle(M.per))
            set(FigProp.h3,'CData',tmp)
        end
    end
end

%% RBC project: extra plots

if Model.Project == 1
    
    %Scale the concentration to the HCT value
    Cblood_inc = Cblood(T_base/Seq.TR+1:end);
    C_scaled = Cblood_inc*Scale_factor;
    Relax2_inc = Relax2(T_base/Seq.TR+1:end);
    Cmax = max(Cblood(:));
    Cmax_scaled = Cmax*Scale_factor;
    
    %In vitro fit, from van Osch et al., 2003 (Figure 4)
    x_in_vitro = [0:0.0001:0.018];
    y_in_vitro = 574451*(x_in_vitro.^2) + 7617.7*x_in_vitro;
    
    figure,
    subplot(2,2,1),
    plot(t,C,'b');ylabel('C (mM)'),xlabel('t (s)'),title('Mean [CA]'),axis square, axis ([0 Model.Tmax 0 Cmax])
    subplot(2,2,2),
    plot(t,Cblood,'b');ylabel('C (mM)'),xlabel('t (s)'),title('[CA] in blood plasma'),axis square, axis ([0 Model.Tmax 0 Cmax])
    subplot(2,2,3),
    plot(t,abs(S),'b');ylabel('|S| (a.u.)'),xlabel('t (s)'),title('Magnitude of signal'),axis square, axis ([0 Model.Tmax 0 inf])
    subplot(2,2,4),
    plot(t,angle(S),'b');ylabel('arg(S)'),xlabel('t (s)'),title('Phase of signal'),axis square, axis ([0 Model.Tmax -inf inf])
    
    figure,
    plot(C_scaled,Relax2_inc,'b','LineWidth',1); hold on; plot(x_in_vitro*1000,y_in_vitro,'r','LineWidth',1), xlabel('C (mM)'), ylabel('R2*(/s)'), title('R2* vs. [CA]'), legend('Simulation','In vitro'), axis square, axis ([0 Cmax_scaled 0 400])
    
    figure,
    plot(C_scaled,Relax2_inc,'b','LineWidth',1); xlabel('C (mM)'), ylabel('R2*(/s)'), title('R2* vs. [CA]'), axis square, axis ([0 Cmax_scaled 0 400])
    
end