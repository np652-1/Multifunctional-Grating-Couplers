function NODES = identical_tree_generation(coordv,coordh,Npv,Nph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            This function generates the binary tree by keep dividing the 
%            domain. The current number of tree level is chosen in such a
%            way that the number of points in the leaf panel is at least
%            8*8 = 64.               
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: 
%      coordv: vertical grid coordinates    
%      coordh: horizontal grid coordinates
% output: NODES sturcture
%     NODES{01,ibox}   the level, staring from 1
%     NODES{02,ibox}   index for four corners [Nhmin,Nhmax,Nvmin,Nvmax]';
%     NODES{03,ibox}   the parent
%     NODES{04,ibox}   the children
%     NODES{05,ibox}   if leaf flag
%     NODES{06,ibox}   x value on the panel, partial list same order as the full x list
%     NODES{07,ibox}   y value on the panel, partial list same order as the full y list
%     NODES{26,ibox}   index of points in the ibox panel in the full Nv*Nh
%                      vector
%     NODES{26,ibox}   index of points in the ibox panel in the full Nv*Nh
%                      vector
%%%%%% misc information
% labeling has nothing to do with nodes
%     NODES{30,1} horizontal grid, vector containing horizontal coordinates
%     NODES{30,2} vertical grid, vector containing vertical coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function is modified from the original tree_generation by assuming
%the leaf panels are identical. As a result, the number of levels are
%automatically determined by the leaf vertical and horizontal panel
%numbers.

% Also, it specializes in metasurface like simulations where Lv << Lh, so
% it generate the tree by making the vertical cut first.

% Npv: vertical leaf panel number
% Nph: horizontal leaf panel number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nv = length(coordv);
    Nh = length(coordh);
    ycutnum = log2(Npv);
    xcutnum = log2(Nph);


nlevels = xcutnum + ycutnum + 1;


 if nlevels <= 1
        error('There is only one level. Make the problem size larger or use direct inversion method instead.')
 end
    % temp variable for generating trees
    BOX = zeros(12,2^(nlevels+1));
    %BOX(1,) = ilevel
    %BOX(2,) = Nhmin
    %BOX(3,) = Nhmax
    %BOX(4,) = Nvmin
    %BOX(5,) = Nvmax
    %BOX(6:9,) = nomeaning
    %BOX(10,) = parent
    %BOX(11,) = son1
    %BOX(12,) = son1


    % determine how many levels I need



    % top level
    BOX(1:12,1) = [1,1,Nh,1,Nv,NaN,NaN,NaN,NaN,NaN,-1,-1]';

    ibox_last = 0;
    ibox_new = 1;

    for ilevel = 1:(nlevels-1)

    ibox_first = ibox_last+1; % starting index for the last level
    ibox_last = ibox_new;   % last index for the last level

        % Loop over all boxes on the level that was last created.
        for ibox = ibox_first:ibox_last

        Lhmin = BOX(2,ibox);
        Lhmax = BOX(3,ibox);
        Lvmin = BOX(4,ibox); 
        Lvmax = BOX(5,ibox);


        if ilevel <= xcutnum % make the vertical cut first
            % make a vertical cut
            Lhhalf = floor((Lhmax + Lhmin)/2);

            box_geom_son1 = [Lhmin,Lhhalf,Lvmin,Lvmax];
            box_geom_son2 = [Lhhalf+1,Lhmax,Lvmin,Lvmax];
        else
            % make a horizontal cut
            Lvhalf = floor((Lvmax + Lvmin)/2);
            box_geom_son1 = [Lhmin,Lhmax,Lvhalf+1,Lvmax];
            box_geom_son2 = [Lhmin,Lhmax,Lvmin,Lvhalf];
        end

        % Create son 1
        ibox_new = ibox_new + 1;
        BOX(1:12,ibox_new) = [ilevel+1,box_geom_son1,NaN,NaN,NaN,NaN,ibox,-1,-1]';
        BOX(11,ibox) = ibox_new;

        % Create son 2
        ibox_new = ibox_new + 1;
        BOX(1:12,ibox_new) = [ilevel+1,box_geom_son2,NaN,NaN,NaN,NaN,ibox,-1,-1]';
        BOX(12,ibox) = ibox_new;
        end

    end


    %%
%     figure;
%     hold on;
%     plot(BOX(2,BOX(1,1:ibox_new)==nlevels) + 1i*BOX(4,BOX(1,1:ibox_new)==nlevels),'x')
%     plot(BOX(3,BOX(1,1:ibox_new)==nlevels) + 1i*BOX(5,BOX(1,1:ibox_new)==nlevels),'x')
%     plot(BOX(2,BOX(1,1:ibox_new)==nlevels) + 1i*BOX(5,BOX(1,1:ibox_new)==nlevels),'x')
%     plot(BOX(3,BOX(1,1:ibox_new)==nlevels) + 1i*BOX(4,BOX(1,1:ibox_new)==nlevels),'x')
%     hold off;

    nboxes  = ibox_new; % number of total boxes over all levels

    NODES   = cell(45,nboxes);
    for ibox = 1:nboxes
        NODES{01,ibox} = BOX(1,ibox);   % the level, staring from 1
        NODES{02,ibox} = BOX(2:5,ibox); % index for four corners
        NODES{03,ibox} = BOX(10,ibox);  % the parent
        NODES{04,ibox} = [];  % the children

        for j = 11:12
            if (BOX(j,ibox) > 0)
                 NODES{04,ibox} = [NODES{04,ibox},BOX(j,ibox)]; % fill in the sons
            end
        end

        NODES{05,ibox} = length(NODES{04,ibox}); % if leaf flag
    end
    
    
    [X,Y] = meshgrid(coordh,coordv);
    XX = X(:);
    YY = Y(:);
    
    for ibox = 1:size(NODES,2)
        if NODES{05,ibox} == 0 % if leaf node
            NNhmin = NODES{02,ibox}(1);
            NNhmax = NODES{02,ibox}(2);
            NNvmin = NODES{02,ibox}(3);
            NNvmax = NODES{02,ibox}(4);

            %leftnum = Nv*(NNhmin-1);
            %%%%% points above
            topnum = (NNvmin-1);
            %%%%% points below
            %belownum = (NNvmax-1)
            %%%%% number each row
            rownum = (NNvmax - NNvmin + 1);
            %%%%% repeat (num of column) times
            colnum = (NNhmax - NNhmin + 1);

            index = [];
            for i = 0:(colnum - 1)
                start = Nv*(NNhmin-1 + i) + topnum;
                I = (start+1):(start + rownum);
                index = [index,I];
                %length((start+1):(start + rownum))
            end
    %%%%%%%%%%%%%

            xvalue2 = XX(index);
            yvalue2 = YY(index);
            %N2 = length(xvalue2);

            NODES{06,ibox} = xvalue2;
            NODES{07,ibox} = yvalue2;
            NODES{26,ibox} = index;
        end
    end
    
    
    % store misc information for convenience
    % labeling has nothing to do with nodes
    %     NODES{30,1} horizontal grid
    %     NODES{30,2} vertical grid
    NODES{30,1} = coordh;
    NODES{30,2} = coordv;
 
end
