function[currBlock] = updateBlock(x,blockData,currBlock,Nblocks)
% Used by getSignedDistanceFunction
% blockData must have Nblocks+1 points, the last two points are the same.
% For points in the last block, blockData(currBlock+1, 1) should exist and
% equal blockData(currBlock, 1)
while(x >= (blockData(currBlock+1,1)-1e-16))
    if (currBlock < Nblocks)
        currBlock = currBlock+1;
    else
        break;
    end
end
end