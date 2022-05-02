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

function FigProp = DisplayFig_gen(G,M,CA,S,C,t)


% Magnetization lattices
figure,
VessContour = contourc(full(double(G.vasc.P)),1);
subplot(2,3,1),
h1 = imagesc(abs(M.per));axis square, axis off,title('abs(Mper)')%;colormap gray
inds = 1;
hold on
while inds < size(VessContour,2)
    subplot(2,3,1),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
    inds = inds + VessContour(2,inds)+1;
end
hold off

subplot(2,3,2),
h2 = imagesc(angle(M.per));axis square, axis off,title('Arg(Mper)')
inds = 1;
hold on,
while inds < size(VessContour,2)
    subplot(2,3,2),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
    inds = inds + VessContour(2,inds)+1;
end
hold off

subplot(2,3,3),
tmp = CA;
tmp(logical(G.vasc.P)) = NaN;
h3 = imagesc(tmp);axis square, axis off,title('C (mM)')
inds = 1;
hold on,
while inds < size(VessContour,2)
    subplot(2,3,3),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
    inds = inds + VessContour(2,inds)+1;
end
hold off

% CA and Signal
subplot(2,3,4),
p1 = plot(0,0);ylabel('C (mM)'),xlabel('t (s)'),title('[CA]'),axis square
subplot(2,3,5),
p2 = plot(0,0);ylabel('|S| (a.u.)'),xlabel('t (s)'),title('Abs(S)'),axis square
subplot(2,3,6),
p3 = plot(0,0);ylabel('arg(S)'),xlabel('t (s)'),title('Arg(S)'),axis square

FigProp.h1 =h1;
FigProp.h2 =h2;
FigProp.h3 =h3;
FigProp.p1 =p1;
FigProp.p2 =p2;
FigProp.p3 =p3;