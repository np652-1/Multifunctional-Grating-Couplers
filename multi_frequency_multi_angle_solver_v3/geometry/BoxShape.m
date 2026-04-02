function[density] = BoxShape(x1,y1,x2,y2,xyGrid,kSharp)
% Create a density distribution of a box shape
BlockData = [x1, y1; x1, y2; x2, y1; x2, y1];
Nblocks = size(BlockData, 1) - 1; % Divide the simulation region along the x direction into Nblock blocks
[distFx] = getSignedDistanceFunction(BlockData, Nblocks, xyGrid);
density = 1 ./ (1 + exp(-2*kSharp*distFx));
end