clear all
close all
clc

gps = importdata('all.gps');
covar = importdata('ukbb.nim.covar');
PC1 = covar.data(:,5);
nimPC = importdata('expandedNIMprojections.txt');
indNa = find(gps(:,1)~=-9);

longitude = gps(indNa,1);
latitude = gps(indNa,2);
nimPC = nimPC(indNa,:);
PC1 = PC1(indNa,:);
n = 20;
qLong = prctile(longitude,[0:100/n:100]);
qLat = prctile(latitude,[0: 100/n:100]);

location = zeros(n,n,2);
pc = zeros(n);
pcAll = zeros(n);
for i = 1:n
    indLong = intersect(find(longitude >= qLong(i)) , find(longitude < qLong(i+1)));
    for j = 1:n
        indLat = intersect(find(latitude >= qLat(j)) , find(latitude < qLat(j+1)));
        ind = intersect(indLong,indLat);
        if length(ind) > 10
            location(i,j,:) =[median(longitude(ind)); median(latitude(ind))];
            pc(i,j) = mean(nimPC(ind,1));
            pcAll(i,j) = mean(PC1(ind));
        end
    end
end

subplot('Position',[0.07, 0.1, 0.45, 0.75])
worldmap([49 59],[-8 2])
geoshow('landareas.shp')
title('PC-1')

ind = find(isnan(pcAll)~=1);
[sorted indsort] = sort(pcAll(ind));
pcAll(ind(indsort)) = [1:length(ind)]/length(ind);
hold on
for i = 1:n
    for j = 1:n
        if isnan(pcAll(i,j))~=1
        circlem(location(i,j,2),location(i,j,1),5,'facecolor',[1-pcAll(i,j) 0.5 1-pcAll(i,j)] )
        end
    end
end

subplot('Position',[0.54, 0.1, 0.45, 0.75])
worldmap([49 59],[-8 2])
geoshow('landareas.shp')
title('NIM PC-1')
hold on
ind = find(isnan(pc)~=1);
[sorted indsort] = sort(pc(ind));
pc(ind(indsort)) = [1:length(ind)]/length(ind);
for i = 1:n
    for j = 1:n
        if isnan(pc(i,j))~=1
        circlem(location(i,j,2),location(i,j,1),5,'facecolor',[1-pc(i,j) 0.5 1-pc(i,j)])
        end
    end
end

subplot('Position',[0.4, 0.9, 0.2, 0.05])
hold off
hold on
for i = 1:100
    plot([i/100 i/100],[ 0.1 0.12], '-','Color', [1-0.01*i 0.5 1-0.01*i],'LineWidth',10);
end
text(-0.1, 0.127, 'Rank (low)','FontSize',16);
text(0.9, 0.127, 'Rank (high)','FontSize',16);

axis off


set(gcf,'PaperPosition',[0 0 20 10])
saveas(1,'gpsNIMPC1.png');

function [h,circlelat,circlelon] = circlem(lat,lon,radius,varargin)
%CIRCLEM draws circles on maps. 
% 
%% Syntax
% 
%  circlem(lat,lon,radius)
%  circlem(...,'units',LengthUnit)
%  circlem(...,'PatchProperty',PatchValue)
%  h = circlem(...)
%  [h,circlelat,circlelon] = circlem(...)
% 
%% Description 
%
% circlem(lat,lon,radius) draws a circle or circles of radius or radii
% given by radius centered at lat, lon, where radius, lat, and
% lon may be any combination of scalars, vectors, or MxN array. All non-
% scalar inputs must have matching dimensions. 
% 
% circlem(...,'units',LengthUnit) specifies a length unit of input
% radius. See validateLengthUnit for valid units. Default unit is
% kilometers. 
% 
% circlem(...,'PatchProperty',PatchValue) specifies patch properties such
% as edgecolor, facecolor, facealpha, linewidth, etc. 
% 
% h = circlem(...) returns the patch handle of plotted circle(s). 
%  
% [h,circlelat,circlelon] = circlem(...) also returns arrays of latitudes 
% and longitudes corresponding to the outline of each circle drawn.  Each
% "circle" is in actuality a polygon made of 100 lat/lon pairs.
% 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas Institute 
% for Geophysics (UTIG) on October 14, 2014.  
%
% See also scircle1, distance, validateLengthUnit
%% Initial input checks: 
assert(license('test','map_toolbox')==1,'circlem requires Matlab''s Mapping Toolbox.')    
try
    MapIsCurrent = ismap; 
    assert(MapIsCurrent==1,'The circlem function requires you to initialize a map first.') 
catch
    error('The circlem function requires you to initialize a map first.') 
end    
assert(nargin>2,'circlem function requires at least three inputs--latitude(s), longitude(s), and radius(ii)')
assert(isnumeric(lat)==1,'Input latitude(s) must be numeric.');
assert(isnumeric(lon)==1,'Input longitude(s) must be numeric.');
assert(max(abs(lat))<=90,'Latitudes cannot exceed 90 degrees North or South.') 
assert(max(abs(lon))<=360,'Longitudes cannot exceed 360 degrees North or South.') 
assert(isnumeric(radius)==1,'Radius must be numeric.') 
%% Declare units
units = 'km'; % kilometers by default
tmp = strncmpi(varargin,'unit',4); 
if any(tmp)
    units = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
end
%% Reshape inputs as needed
% columnate:
lat = lat(:); 
lon = lon(:); 
radius = radius(:); 
% How many circles do we need to make? 
NumCircles = max([length(lat) length(lon) length(radius)]); 
% Make vectors all the same lengths so scircle1 will be happy:  
if length(lat)<NumCircles
    assert(isscalar(lat)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    lat = lat*ones(NumCircles,1); 
end
if length(lon)<NumCircles
    assert(isscalar(lon)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    lon = lon*ones(NumCircles,1); 
end
if length(radius)<NumCircles
    assert(isscalar(radius)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    radius = radius*ones(NumCircles,1); 
end
%% Calculate circle coordinates:
[circlelat,circlelon] = scircle1(lat,lon,radius,[],earthRadius(units));
%%  Plot and format circle(s): 
h = patchm(circlelat,circlelon,'k','facecolor','none'); 
if nargin>3 && ~isempty(varargin)
    set(h,varargin{:})
end
%% Clean up: 
if nargout==0
    clear h
end
end