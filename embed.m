function [ coordinatesMatrix ] = embed( adj,o,l,b,x,t,r,u,e,L,i)
%EMBED Summary of this function goes here
%  Inputs:
%           adj: adjacency matrix of the graph
%           o: path for the output file
%           l: path for the Dist file
%           b: number of landmarks that make a first calibration of their
%           coordinates                                   
%           x: the number of dimensions
%           t: input file path
%           r: path for the landmarks file
%           e: curvature of hyperbolic space
%           i: number of landmarks with which each node calibrates its
%           coordinates
% Outputs:                                                                 
%           oordinatesMatrix: A n x m matrix in which every row represents
%           a node and every column the coordinate of the node in that
%           dimension (n nodes, m dimensions total)
% Author:   Konstantinos Tsitseklis.


% 'inputs\scF\scF100' is the default file destination. For the code to work either you have to create
% these directories or change it with one of your choosing.

  landmarkSelect(adj,'',str2num(L),'.inputs\scF\scF100');
  writeFiles(o,l,b,x,t,r,u,e,L,i)
  
    !rigel2.exe -1 
    !rigel2.exe 0  
    
    land = importdata('outputs\scF\scF100.land');
    nodes = importdata('outputs\scF\scF1000.coord');
    
    coordinatesMatrix = [ land; nodes];
    coordinatesMatrix = sortrows( coordinatesMatrix );
end

