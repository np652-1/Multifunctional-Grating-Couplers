function NODES = near_surround_proxy_surface_creation(NODES,proxy_layer,Nph)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           This function construct proxy surface made of points surrounding a panel
%           empty space represents the panel'o' represents points of the proxy
%           surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           here is an illustration an one-layer-thickness proxy surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    top line
%                                  ooooooooooo
%                                  o         o
%                                  o         o
%                             left o  panel  o right
%                             line o  points o line
%                                  o  are    o
%                                  o  inside o
%                                  o         o
%                                  o         o
%                                  ooooooooooo
%                                   bottom line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Input:
%                NODES structure, with binary tree created (tree_generation)
%                                 and x,y position coordinates stored in
%                                 each panel.              (tree_generation)
%                proxy_layer: the number of layer surrounding each panel.
%
%                According to Table 1 in the Gopal Martinsson paper,
%                choose proxy_layer = 1 for compression accuracy > 1e-5
%                choose proxy_layer = 2 for compression accuracy > 1e-10
%                choose proxy_layer = 3 for compression accuracy > mach_eps
%                
%
%
%           Outputs:
%                  NODES sturcture, with x,y position coordinates for the proxy
%                            surface for each panel computed and stored
%
%                          NODES{13,ibox}   x coordinates of proxy surface
%                          NODES{14,ibox}   y coordinates of proxy surface
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nbox = size(NODES,2);
    xcutnum = log2(Nph);
    Nh = length(NODES{30,1}); % horizontal grid number
    Nv = length(NODES{30,2}); % vertical grid number
    h = NODES{30,1}(2) - NODES{30,1}(1);
    c = 5; % distance paramter of the far away proxy surface
    outnum = c*Nv;
    for ibox = 2: nbox % exclude root node
        level = NODES{1,ibox};
        
    if ibox == 2^(level-1) %first of the level
        % orientation of horizontal axis
        if (NODES{30,1}(end) - NODES{30,1}(1)) > 0
            Oh = 1; %positive orientation
        elseif ((NODES{30,1}(end) - NODES{30,1}(1)) < 0)
            Oh = -1; %negative orientation
        else 
            Oh = 0; % will return error
        end


        % orientation of vertical axis
        if (NODES{30,2}(end) - NODES{30,2}(1)) > 0
            Ov = 1; %positive orientation
        elseif (NODES{30,2}(end) - NODES{30,2}(1)) < 0
            Ov = -1; %negative orientation
        else 
            Ov = 0; % will return error
        end


        %Oh = (NODES{30,1}(end) - NODES{30,1}(1)) > 0;
        %Ov = (NODES{30,2}(end) - NODES{30,2}(1)) > 0;

        % index to position conversion is then determined by orientation


        Nhmin = NODES{02,ibox}(1);
        Nhmax = NODES{02,ibox}(2);
        Nvmin = NODES{02,ibox}(3);
        Nvmax = NODES{02,ibox}(4);

        
        if level > (xcutnum + 1) % surrounding proxy surface
        
        % horizontal line x index
        indxh = (Nhmin-proxy_layer):(Nhmax+proxy_layer);
        lenx = length(indxh);
        indxh = repmat(indxh,proxy_layer,1);
        % vertical line y index
        indyv = Nvmin:Nvmax;
        leny = length(indyv);
        indyv = repmat(indyv,proxy_layer,1);


        % horizontal line y index
            % line above
        indyh_up = (Nvmax+1:(Nvmax+1 + proxy_layer - 1))';
        indyh_up = repmat(indyh_up,1,lenx); 
            % line below
        indyh_down = ((Nvmin-1 -proxy_layer +1):Nvmin-1)';
        indyh_down = repmat(indyh_down,1,lenx);


        % vertical line x index
            % left line
        indxv_left = (Nhmin-1 - proxy_layer +1):Nhmin-1;
        indxv_left = repmat(indxv_left,leny,1)';
            % right line
        indxv_right = Nhmax+1:(Nhmax+1 + proxy_layer  - 1);
        indxv_right = repmat(indxv_right,leny,1)';

        %[top;down;left;right]
        indxh = indxh(:);
        indxv_left = indxv_left(:);
        indxv_right = indxv_right(:);


        if Oh == 1 % list of hori coordinates is in increasing order
            xp = ([indxh;indxh;indxv_left;indxv_right]-1)*h;
        elseif Oh == -1 % list of hori coordinates is in decreasing order
            xp = (Nh-[indxh;indxh;indxv_left;indxv_right])*h;
        else
            error('Horizontal orientation invalid.');
        end

        indyv = indyv(:);
        indyh_up = indyh_up(:);
        indyh_down = indyh_down(:);



        if Ov == 1 % list of vertical coordinates is in increasing order
            yp = ([indyh_up;indyh_down;indyv;indyv]-1)*h;
        elseif Ov == -1 % list of vertical coordinates is in decreasing order
            yp = (Nv-([indyh_up;indyh_down;indyv;indyv]))*(h);
        else
            error('Horizontal orientation invalid.');
        end


        NODES{13,ibox} = xp;  %x coordinates of proxy surface
        NODES{14,ibox} = yp;  %y coordinates of proxy surface


        % figure;
        % hold on;
        % plot(NODES{06,ibox}+1i*NODES{07,ibox},'o');
        % % proxy points
        % plot(xp+1i*yp,'*');
        % hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif  level <= (xcutnum + 1) % vertical line sticking out proxy surface
        % horizontal line x index
        indxh = (Nhmin-proxy_layer - outnum):(Nhmax+proxy_layer + outnum);
        lenx = length(indxh);
        indxh = repmat(indxh,proxy_layer,1);
        
        % vertical line y index
        indyv = (Nvmin-outnum):(Nvmax+outnum);
        leny = length(indyv);
        indyv = repmat(indyv,proxy_layer,1);
        

        % horizontal line y index
            % line above
        indyh_up = ((Nvmax+1+outnum):(Nvmax+1 + proxy_layer - 1+outnum))';
        indyh_up = repmat(indyh_up,1,lenx); 
            % line below
        indyh_down = ((Nvmin-1 -proxy_layer +1 - outnum):(Nvmin-1-outnum))';
        indyh_down = repmat(indyh_down,1,lenx);


        % vertical line x index
            % left line
        indxv_left = (Nhmin-1 - proxy_layer +1 - outnum):(Nhmin-1-outnum);
        indxv_left = repmat(indxv_left,leny,1)';
            % right line
        indxv_right = (Nhmax+1+outnum):(Nhmax+1 + proxy_layer  - 1+outnum);
        indxv_right = repmat(indxv_right,leny,1)';

        %[top;down;left;right]
        indxh = indxh(:);
        indxv_left = indxv_left(:);
        indxv_right = indxv_right(:);
        
        proxy_layer_vert = 4;
%%%%%%% near coordinates
        % vertical line y index near
        indyvn = Nvmin:Nvmax;
        lenyn = length(indyvn);
        indyvn = repmat(indyvn,proxy_layer_vert,1);
        
        % vertical line x index
            % left line
        indxvn_left = (Nhmin-1 - proxy_layer_vert +1):Nhmin-1;
        indxvn_left = repmat(indxvn_left,lenyn,1)';
            % right line
        indxvn_right = Nhmax+1:(Nhmax+1 + proxy_layer_vert  - 1);
        indxvn_right = repmat(indxvn_right,lenyn,1)';        
        
        
        indxvn_left = indxvn_left(:);
        indxvn_right = indxvn_right(:);
        
        
        indyvn = indyvn(:);

        if Oh == 1 % list of hori coordinates is in increasing order
            xp = ([indxh;indxh;indxv_left;indxv_right;indxvn_left;indxvn_right]-1)*h;
        elseif Oh == -1 % list of hori coordinates is in decreasing order
            xp = (Nh-[indxh;indxh;indxv_left;indxv_right;indxvn_left;indxvn_right])*h;
        else
            error('Horizontal orientation invalid.');
        end
        
        indyv = indyv(:);
        indyh_up = indyh_up(:);
        indyh_down = indyh_down(:);





        if Ov == 1 % list of vertical coordinates is in increasing order
            yp = ([indyh_up;indyh_down;indyv;indyv;indyvn;indyvn]-1)*h;
        elseif Ov == -1 % list of vertical coordinates is in decreasing order
            yp = (Nv-([indyh_up;indyh_down;indyv;indyv;indyvn;indyvn]))*(h);
        else
            error('Horizontal orientation invalid.');
        end


        NODES{13,ibox} = xp;  %x coordinates of proxy surface
        NODES{14,ibox} = yp;  %y coordinates of proxy surface
        end
        
    end    
    end

end