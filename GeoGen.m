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

function P=GeoGen(Model)

% Generate the geometry matrices
% Geometry are stored/load in Model.geo.path directory if any.
%
% Input:
%   - geo: structure containing geometry parameters
%
% Output:
%   - P: structure containing geometry matrices

% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012


%% Initialize
vasc.P  = false( Model.geo.res );          %matrice de position
vasc.Pp = false( Model.geo.res );          %matrice de position
cell.P  = false( Model.geo.res );
if Model.Project == 1
    fov = 70e-6;
else
    fov = Model.geo.vasc.R * sqrt( Model.geo.vasc.N * pi / Model.geo.vasc.vp );
end
verbose = Model.Flag.verbose;

%% Generate Geo filename
vasc.filename = sprintf('GeoV_v%d_res%d_R%d_N%d_Int%.1f_Vp%d_Ww%d',Model.geo.vasc.Id,Model.geo.res,Model.geo.vasc.R*1e6, ...
    Model.geo.vasc.N,Model.geo.vasc.Inter*1e6,Model.geo.vasc.vp*1e2,Model.geo.vasc.wallw*1e6,Model.geo.vasc.Id);
cell.filename = sprintf('GeoC_v%d_c%d_res%d_R%d_N%d_Int%.1f_Vp%d_Ww%d_P%d_Rc%d_CInter%.f',Model.geo.vasc.Id, ...
    Model.geo.cell.Id,Model.geo.res,Model.geo.vasc.R*1e6,Model.geo.vasc.N,Model.geo.vasc.Inter*1e6,Model.geo.vasc.vp*1e2,Model.geo.vasc.wallw*1e6,Model.geo.cell.Poro*1e2,Model.geo.cell.Rc*1e6,Model.geo.cell.Inter*1e6);
ind = regexp(vasc.filename,'[.]');
vasc.filename(ind) = 'p';
vasc.filename = fullfile(Model.geo.path,vasc.filename);
vasc.filename = [vasc.filename '.mat'];
ind = regexp(cell.filename,'[.]');
cell.filename(ind) = 'p';
cell.filename = fullfile(Model.geo.path,cell.filename);
cell.filename = [cell.filename '.mat'];

%% Geo already generated ?
if(exist(vasc.filename,'file')~=0)        % for the vessels ?
    if Model.Flag.verbose,disp('Loading vessel geometry...');end
    load(vasc.filename);
    P.vasc = vasc;
    if(exist(cell.filename,'file')~=0)    % for the cells ?
        if Model.Flag.verbose,disp('Loading cell geometry...');end
        load(cell.filename);
        P.cell = cell;
    else  % Generate the cells                               
        if Model.Flag.verbose,disp('Generating cell geometry...');end
        model = CellGenerator(Model.geo.cell.Rc*1e6,Model.geo.cell.Inter*1e6,Model.geo.cell.Poro,vasc.radii*1e6,vasc.center*1e6,fov*1e6,vasc.indices,verbose);
        cell.radii = model.Radius*1e-6;
        cell.center = model.Centers*1e-6;
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%% Parametres geometriques
%        dspatial = fov ./ Model.geo.res;                   % vector with size of the subvoxels. 
        dx=fov/Model.geo.res;                   % Resolution spatiale en y
        dy=fov/Model.geo.res;                   % Resolution spatiale en y

        %Echantillonage spatiale
        [x,y] = ndgrid( 1:Model.geo.res, 1:Model.geo.res );   %Creation de trois vecteur n*n*n rempli avec des valeurs allant de -n/2 a n/2-1 en variant de 1 seulement dans la direction du vecteur
        %x = dspatial * x;
        x = dx * x;
        y = dy * y;
        xy = [x(:)  y(:)]';  
        Tmp2 = zeros(Model.geo.res,Model.geo.res);
        for NbEl = 1:numel(cell.radii)
            Tmp = zeros(Model.geo.res,Model.geo.res);
            a = cell.center(NbEl,1);                %Position initiale du centre des cylindres
            b = cell.center(NbEl,2);
            R = cell.radii(NbEl);
            
            r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);
            indx_inside = r<R;                        %indice des pixels a l'interieur du cylindre
            Tmp(indx_inside)= 1;
            
            Tmp2 = Tmp + Tmp2;
        end
        
        clear x y r indx_inside
        
        cell.P = sparse(Tmp2-vasc.P);
        if Model.Flag.SaveGeo, save(cell.filename,'cell');end
        P.cell = cell;
    end
    P.ees.P = sparse(ones(size(P.vasc.P)) - full(P.vasc.P + P.cell.P));
else
    if Model.Flag.verbose,disp('Generating vessel geometry...');end
    if Model.Project == 1
        rand('twister', sum(100*clock)); %sum(0) for fixed positioning of RBCs
    else
        rand('twister', sum(100*clock)); %permet a rand de demarrer chaque fois avec une valeur differente
    end
    if(Model.geo.vasc.Inter ~= -1)
        model = VesselGenerator(Model.Project,Model.Shape,Model.geo.vasc.vp,Model.geo.vasc.R*1e6,Model.geo.vasc.R2*1e6,Model.geo.vasc.N,Model.geo.vasc.Inter*1e6,verbose);
        vasc.radii = model.Radius*1e-6;
        vasc.center = model.Centers*1e-6;
        vasc.indices = model.Indices;
    elseif(mod(Model.geo.vasc.N^(1/2),1) == 0) % vessels on a cartesian grid
        radii = ones([1 Model.geo.vasc.N])*Model.geo.vasc.R;
        tmp = (fov/sqrt(Model.geo.vasc.N)/2:fov/sqrt(Model.geo.vasc.N):fov);
        center(:,1) = tmp;
        center = repmat(center,[sqrt(Model.geo.vasc.N) 1]);
        for aaa=1:sqrt(Model.geo.vasc.N)
            center((aaa-1)*sqrt(Model.geo.vasc.N)+1:aaa*sqrt(Model.geo.vasc.N),2) = ones([1 sqrt(Model.geo.vasc.N)])*tmp(aaa);
        end
%         model.Radius = radii*1e-6;
%         model.Centers = center*1e-6;
%         model.CubeSize(1) = x_fov*1e-6;
        vasc.radii = radii;
        vasc.center = center;
        vasc.indices = (1:Model.geo.vasc.N);
    else
        if Model.Flag.verbose,disp('Wrong Rext or Nb cylinder values'),end
        return;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Parametres geometriques
    dx=fov/Model.geo.res;                   % Resolution spatiale en x
    dy=fov/Model.geo.res;                   % Resolution spatiale en y

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% Creation matrice geometrie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Pxy=zeros(x_pix,y_pix);          %matrice de position
%     PxyCell = zeros(x_pix,y_pix);

    %Echantillonage spatiale
    [x,y] = ndgrid(1:Model.geo.res,1:Model.geo.res);   %Creation de trois vecteur n*n*n rempli avec des valeurs allant de -n/2 a n/2-1 en variant de 1 seulement dans la direction du vecteur
    x = dx * x;
    y = dy * y;
    xy = [x(:)  y(:)]';                          %Creation vecteur 3*(n^3)

    for it=1:size(vasc.radii,2)
        a = vasc.center(it,1);
        b = vasc.center(it,2);
        Ray = vasc.radii(it);
        if Model.Project == 1
            if Model.Shape == 1     %Ellipsoidal RBCs
                Rax = Model.geo.vasc.R2;
                ry = xy(2,:) - b;
                rx = xy (1,:) - a;
                r_ellipse = sqrt(ry.^2./Ray.^2+rx.^2./Rax.^2);
                indx_RBC = r_ellipse<1;
                vasc.P(indx_RBC)= 1;
            else
                r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);     %Circular RBCs
                indx_inside = r<Ray;
                vasc.P(indx_inside)= 1;
            end
        else
            r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);         %Vessels
            indx_inside = r<Ray;
            vasc.P(indx_inside)= 1;
        end
    end
    
    if Model.Flag.verbose,disp('Generating cell geometry...');end
    model = CellGenerator(Model.geo.cell.Rc*1e6,Model.geo.cell.Inter*1e6,Model.geo.cell.Poro,vasc.radii*1e6,vasc.center*1e6,fov*1e6,vasc.indices,verbose);
    cell.radii = model.Radius*1e-6;
    cell.center = model.Centers*1e-6;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%% Parametres geometriques
    dx=fov/Model.geo.res;                   % Resolution spatiale en x
    dy=fov/Model.geo.res;                   % Resolution spatiale en y
    
    %Echantillonage spatiale
    [x,y] = ndgrid(1:Model.geo.res,1:Model.geo.res);   %Creation de trois vecteur n*n*n rempli avec des valeurs allant de -n/2 a n/2-1 en variant de 1 seulement dans la direction du vecteur
    x = dx * x;
    y = dy * y;
    xy = [x(:)  y(:)]';
    Tmp2 = zeros(Model.geo.res,Model.geo.res);
    for NbEl = 1:numel(cell.radii)
        Tmp = zeros(Model.geo.res,Model.geo.res);
        a = cell.center(NbEl,1);                %Position initiale du centre des cylindres
        b = cell.center(NbEl,2);
        R = cell.radii(NbEl);
        
        r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);
        indx_inside = r<R;                        %indice des pixels a l'interieur du cylindre
        Tmp(indx_inside)= 1;
        
        Tmp2 = Tmp + Tmp2;
    end
    
    clear x y r indx_inside
    
    %Different EES mask for RBC project
    if Model.Project == 1
        ees.P = ones(size(vasc.P)) - full(vasc.P + cell.P);
        cell.P = zeros(size(vasc.P));
    else
        cell.P = (Tmp2-vasc.P);
    end
    
    vasc.P = (vasc.P);

    %% %%%%%%%%%%%%%%%%%% selection couronne %%%%%%%%%%%%%%%%%%%%%
    Tmp = conv2(full(vasc.P),[0 1 0; 1 1 1; 0 1 0],'same');
    Tmp(Tmp > 1) = 1;
    vasc.Pp=Tmp-full(vasc.P);             %couronne de 1 pixel de 1 d'epaisseur autour du cylindre avec tout le reste a 0
    vasc.Pp = sparse(vasc.Pp);
    
    if Model.Flag.SaveGeo,save(vasc.filename,'vasc');end
    if Model.Flag.SaveGeo,save(cell.filename,'cell');end
        
    % Build geometry structure in sparse matrix for saving
    P.vasc = vasc;
    P.cell = cell;
    if Model.Project == 1
        P.ees.P = ees.P;
    else
        P.ees.P = sparse(ones(size(P.vasc.P)) - full(P.vasc.P + P.cell.P));
    end

end

