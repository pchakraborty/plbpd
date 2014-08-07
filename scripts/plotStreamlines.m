function [] = plotStreamlines(infile, Y, startx, startz)
%function [] = plotStreamlines(hdf5datafile, Y, startx, startz)
%
%  function to plot streamlines in a given y=Y plane
%  the velocity is taken from the input HDF5 file 'infile'
%  startx, startz define the starting positions of the streamlines

if nargin < 4
    disp('ERROR: incorrect number of arguments')
    return
end

% read data
u = hdf5read(infile,'u');
nrows = size(u,2);
n2    = size(u,3); 
ncols = size(u,4);
assert((Y>0)&&(Y<=n2));

[x,z] = meshgrid(1:ncols,1:nrows);
ux = zeros(nrows,ncols);
uz = zeros(nrows,ncols);
speed = zeros(nrows,ncols);
for i=1:nrows
    for j=1:ncols
        ux(i,j) = u(1,j,Y,i);
        uz(i,j) = u(3,j,Y,i);
        speed(i,j) = sqrt(ux(i,j)^2 + uz(i,j)^2);
    end
end

close all, figure(1)
h = streamline(x,z,ux,uz,startx,startz);

