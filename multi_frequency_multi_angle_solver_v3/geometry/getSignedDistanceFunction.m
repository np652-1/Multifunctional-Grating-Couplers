function[distFx,blockData,simData,inMat] = getSignedDistanceFunction(blockData,Nblocks,xySim)

Nsim = size(xySim,1);
inMat = zeros(Nsim,1);
distFx = zeros(Nsim,1);
simData = zeros(Nsim,2);

currBlock = 1;
for simPt=1:Nsim
    % find which block my simPt is in (assumes the simPos x-values are non-decreasing)
    currBlock = updateBlock(xySim(simPt,1),blockData,currBlock,Nblocks);
    simData(simPt,:) = [simPt currBlock]; 
    
    % the current block is made of up to three line segments; get the
    % distances from each
    dist1 = Inf; dist2 = Inf; dist3 = Inf; dist4 = Inf;
    if (currBlock > 1) % there is a "left wall"
        xLine1 = [blockData(currBlock,1); blockData(currBlock,1)];
        yLine1 = sort(blockData(currBlock-1:currBlock,2));
        dist1 = getDistVertLineSegment(xySim(simPt,:),xLine1,yLine1);
    end
    if (currBlock < Nblocks) % there is a "right wall"
        xLine2 = [blockData(currBlock+1,1); blockData(currBlock+1,1)];
        yLine2 = sort(blockData(currBlock:currBlock+1,2));
        dist2 = getDistVertLineSegment(xySim(simPt,:),xLine2,yLine2);
    end
    % there is always a horizontal wall, whose distance is just abs(height-y)
    if currBlock>1 && currBlock<Nblocks
        dist3 = abs(blockData(currBlock,2) - xySim(simPt,2));
        dist4 = abs(blockData(1,2) - xySim(simPt,2));
    end
    
    % inMat = +1 if height > y, else -1
    if currBlock>1 && currBlock<Nblocks
        inMat(simPt) = sign(blockData(currBlock,2) - xySim(simPt,2)) * sign(xySim(simPt,2) - blockData(1,2));
    else
        inMat(simPt) = -1;
    end
    distFx(simPt) = inMat(simPt) * min([dist1,dist2,dist3,dist4]);
end
end