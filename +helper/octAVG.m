function [RESPavg, BFavg, lHan] = octAVG(BF,RESP,MINpts,Step_oct, plotVar, Wind_oct, BFmin)
% File octAVG.m
%
% Calculates a smoothed curve from POPULATION DATA
% Smoothed curve is a triangular-weighted AVG over 'Wind_oct' octaves
% A point is calculated every 'Step_oct' octaves
%

if ~exist('MINpts','var')
   MINpts=2;   % Default: require at least 2 points in window
end

if ~exist('Step_oct','var')
   Step_oct=.25;   % Default: 1/2-octave smoothing, 1/4-octaves steps
end
if ~exist('plotVar', 'var')
    plotVar= 0;
end

BF= BF(:);
RESP= RESP(:);

if ~exist('Wind_oct', 'var')
    Wind_oct=2*Step_oct;
elseif isempty(Wind_oct)
    Wind_oct=2*Step_oct;
end

if max(BF)<20
    CFscale= 1e3;
else 
    CFscale= 1;
end
if ~exist('BFmin', 'var')
    BFmin=100/CFscale;
end
BFmax=20e3/CFscale;

BFstep=BFmin;
BFind=1;
while(BFstep<=BFmax)
   BFavg(BFind)=BFstep;
   Winds=find((BF>=BFstep/(2^(Wind_oct/2)))&(BF<=BFstep*(2^(Wind_oct/2)))&(~isnan(RESP)));  % Remove any NaNs in RESP
   if length(Winds)>=MINpts
      weights=1-abs(log2(BF(Winds)/BFstep))/(Wind_oct/2); % triangular window
      weights=weights.*(weights>0);  % Set all weights outside window to 0 (i.e., rectify for this definition)
      RESPavg(BFind)=1/sum(weights)*sum(weights.*RESP(Winds));
   else
      RESPavg(BFind)=NaN;
   end
   %% Move to next BFstep
   BFstep=BFstep*2^Step_oct;
   BFind=BFind+1;
end

if plotVar
    lHan= plot(BFavg, RESPavg, '-', 'linew', 2);
else
    lHan= nan;
end
    