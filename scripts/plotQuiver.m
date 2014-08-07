#!/usr/bin/octave -qf
% function [] = plotQuiver(filename, Lx, Lz, scale)
function [] = plotQuiver(filename, Lx, Lz, scale)

velocity = load(filename);
ux = velocity(1:Lz,1:Lx);
uz = velocity(2*Lz+1:3*Lz,1:Lx);
[x, z] = meshgrid(1:1:Lx,1:1:Lz);
quiver(x,z,ux,uz,scale);
