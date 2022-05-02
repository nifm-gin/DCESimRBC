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

function Cp = AIFGen(AIF,t,dt)

% Return the value of the vascular CA concentration
% If t is an array, Cp is an array.
%
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012


switch lower(AIF.mode)
    case {'slow'} % slow bolus (Beaumont 2009)
        if numel(AIF.Ta) > 1 || numel(AIF.Cpic)>1
            disp('Cpic and/or Ta is badly conditionned (need only 1 element for this aif mode)');
        else
            Tp=AIF.Ta+30;           % rise time (s)
            T0 = AIF.Ta;           % arrival time (s)
            Cpic = AIF.Cpic;       % concentration at peak
            if any(t < T0), Cp(t < T0) = 0; end
            if any(t >= T0 & t <= Tp), Cp(t >= T0 & t <= Tp) = Cpic/(Tp-T0) * (t(t >= T0 & t <= Tp) - T0); end
            if any(t > Tp),Cp(t > Tp) = Cpic*(1.94/2.7*exp(-(t(t > Tp)-Tp)*(1.01/60)) + 0.76/2.7*exp(-(t(t > Tp)-Tp)*(0.03/60)));end
        end
        
    case {'fast'} % first passage bolus
        if numel(AIF.Ta) > 1 || numel(AIF.Cpic)>1
            disp('Cpic and/or Ta is badly conditionned (need only 1 element for this aif mode)');
        else
            K=7000;             % Gamma variate parameter
            T0=AIF.Ta;      % " (s)
            alpha=2.6;          % "
            beta=0.56;          % "
            Tdilu=5;            % dilution time (s)
            DiluAmp=0.2;        % (%) percentage of Cpic
            para = [K T0 alpha beta Tdilu DiluAmp];
            if any(t < T0),Cp(t < T0) =0; end
            if any(t >= T0),Cp(t >= T0) = gammaV_dilu(para,t(t >= T0))/gammaV_dilu(para,alpha*beta+T0)*AIF.Cpic;end
        end
        
    case {'mix'}
        if numel(AIF.Ta) > 1 || numel(AIF.Cpic)>1
            disp('Cpic and/or Ta is badly conditionned (need only 1 element for this aif mode)');
        else
            K=7000;             % Gamma variate parameter
            T0=AIF.Ta;      % " (s)
            alpha=6.6;          % "
            beta=0.26;          % "
            Tdilu=5;            % dilution time (s)
            DiluAmp=0.2;        % (%) percentage of Cpic
            para = [K T0 alpha beta Tdilu DiluAmp];
            if any(t < T0), Cp(t < T0) =0; end
            if any(t >= T0), Cp(t >= T0) = gammaV_mix(para,t(t >= T0))/gammaV_mix(para,alpha*beta+T0)*AIF.Cpic; end
        end
        
    case {'step'}
        if numel(AIF.Ta) > 1 || numel(AIF.Cpic)>1
            disp('Cpic and/or Ta is badly conditionned (need only 1 element for this aif mode)');
        else
            T0 = AIF.Ta;
            if any(t < T0), Cp(t < T0) = 0; end
            if any(t >= T0),Cp(t >= T0) = AIF.Cpic;end
        end
        
    case {'dirac'}
        if numel(AIF.Ta) > 1 || numel(AIF.Cpic)>1
            disp('Cpic and/or Ta is badly conditionned (need only 1 element for this aif mode)');
        else
            T0 = AIF.Ta;
            if any(t<T0), Cp(t <T0) = 0; end
            if any(t >= T0 & t < (T0 + dt)),Cp(t >= T0 & t < (T0 + dt)) = AIF.Cpic; end
            if any(t >= T0+dt), Cp(t >= T0+dt) = 0;end
        end
        
    case {'ladder'}
        if numel(AIF.Ta) ~= numel(AIF.Cpic)
            disp('mode ''Ladder'' is badly conditionned');
        else
            Ta = AIF.Ta;
            if any(t<AIF.Ta(1)), Cp(t<AIF.Ta(1)) = 0; end
            for a=1:numel(AIF.Ta)-1
                Cp(Ta(a) < t & t <= Ta(a+1)) = AIF.Cpic(a);
            end
            Cp(t >= Ta(end)) = AIF.Cpic(end);
        end
        
    otherwise
        disp('Unknown AIF mode');
        Cp = 0;
        return
end
end


function out = gammaV_dilu(para,t)
K=para(1);
T0=para(2);
alpha=para(3);
beta=para(4);
Tdilu=para(5);
DiluAmp=para(6);

out = K*(t - T0).^alpha .* exp( - (t - T0)/beta) + gammaV(para,alpha*beta+T0)*DiluAmp*(1 - exp(-(t - T0)/Tdilu));
end

function out = gammaV(para,t)
K=para(1);
T0=para(2);
alpha=para(3);
beta=para(4);
Tdilu=para(5);
baseline=para(6);

out = K*(t - T0).^alpha .* exp( - (t - T0)/beta);
end

function out = gammaV_mix(para,t)
K=para(1);
T0=para(2);
alpha=para(3);
beta=para(4);
Tdilu=para(5);
DiluAmp=para(6);

out = (K*(t - T0).^alpha .* exp( - (t - T0)/beta) + gammaV(para,alpha*beta+T0)*DiluAmp*(1 - exp(-(t - T0)/Tdilu))).*(1.94/2.7*exp(-t*(1.01/60)) + 0.76/2.7*exp(-t*(0.03/60)));
end