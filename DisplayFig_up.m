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

function FigProp = DisplayFig_up(G,M,CA,S,C,t,FigProp,tt)


set(FigProp.h1,'CData',abs(M.per))
set(FigProp.h2,'CData',angle(M.per))
tmp = CA;
tmp(logical(G.vasc.P)) = NaN;
set(FigProp.h3,'CData',tmp)
if ~isempty(FigProp.p1) && numel(t) > 0
    set(FigProp.p1,'XData',t);
    set(FigProp.p1,'YData',C);
    set(FigProp.p2,'XData',t);
    set(FigProp.p2,'YData',abs(S));
    set(FigProp.p3,'XData',t);
    set(FigProp.p3,'YData',angle(S));
    
    set(get(FigProp.p1,'Parent'),'XLim',[0 tt])
    set(get(FigProp.p2,'Parent'),'XLim',[0 tt])
    set(get(FigProp.p3,'Parent'),'XLim',[0 tt])
else
    set(get(FigProp.p1,'Parent'),'XLim',[0 tt])
    set(get(FigProp.p2,'Parent'),'XLim',[0 tt])
    set(get(FigProp.p3,'Parent'),'XLim',[0 tt])
    
end
drawnow