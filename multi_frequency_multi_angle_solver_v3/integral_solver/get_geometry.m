function NODES = get_geometry(NODES,BB)
% store geometry information in the NODES structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        loop over leaf nodes, get and store B = -k^2*Chi on each panel, which
%        contains the scatterer information.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Input:
%                NODES structure, with binary tree created (tree_generation)
%                                 and x,y position coordinates stored in
%                                 each panel.              (tree_generation)
%
%
%                BB: BB = B(:) is the vertorized version of B = -k^2*Chi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Input:
%                NODES structure, with BB value stored on each leaf panel
%
%                               NODES{10,ibox}      BB value on each leaf panel
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nlevels = NODES{1,end};
    % loop over leaf nodes
    for ibox = 2^(nlevels-1):(2^nlevels - 1)
    %for ibox = 1:size(NODES,2)
        %if NODES{05,ibox} == 0 % if leaf node
            index = NODES{26,ibox};
            n = length(index);
            NODES{10,ibox} = spdiags(BB(index),0,n,n);
            
        %end
    end
end

