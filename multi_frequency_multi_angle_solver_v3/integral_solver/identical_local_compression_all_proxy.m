function NODES = identical_local_compression_all_proxy(NODES,k,h, param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           This function computes the interpolative decomposition (ID) of
%           the off-diagonal Green's function hierarchically.
%           It does so by first using ID to get skeleton points on the leaf nodes,
%           then combining skeleton points of children as points in the
%           parental nodes and compressing those points again by ID to get
%           further compressed skeleton points.
%
%           Proxy surfaces are always used for the lowest level. Then it
%           decides automatically on higher levels whether to use proxy
%           surface, depending on the size of total skeleton points and
%           proxy surfaces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Input:
%                NODES structure, with binary tree created (tree_generation)
%                                 and x,y position coordinates stored in
%                                 each panel,              (tree_generation)
%                                 and location of proxy surface point
%                                 coordinates stored.       (proxy_surface_creation)
%                k: wave number
%                h: grid spacing
%                param (direct solver parameters):
%                       param.compress_acc: relative accuracy for ID
%                       param.print_flag: whether to print compressed rank info, 1 yes, 0 no
%                       param.fileID: where printouts go
%                       param.print_flag: whether to plot skeleton points,1 yes, 0 no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Output:
%                 NODES structure, with compressed off-diagonal Greens's function stored
%                                NODES{08,ibox}   row ID index: J
%                                NODES{09,ibox}   row ID matrix: U
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             off diagonal compression form:
%             |U*Goff12(J1,J2)*(V.') -    Goff12|/|Goff12| < param.compress_acc in Frobenius norm.
%
%             Because of reciprocity, U = V, and J1 = J2, that's why we only need
%             to compute the row ID, not the column ID.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function is modified from the original local_compression by
%  assuming the leaf panels are identical. As a result it ONLY compress the
%  first panel on each level, then copy the compression information to
%  others on that level.
%
%  This function uses proxy surfaces for the compression of all levels.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlevels = NODES{1,end}; % number of levels

% loop over all level
for level = nlevels:-1:1
    %for level = nlevels:-1:(nlevels-1)
    %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
    if param.plot_flag == 1
        figure;
        hold on;
    end
    % loop over all nodes on a level
    %ktot = 0;
    for iibox = 2^(level-1):(2^level-1)
        
        if NODES{05,iibox} == 0 % if it is a leaf node
            
            if iibox == 2^(level-1) % only compress the first elements
                xvalue_rest = NODES{13,iibox};
                yvalue_rest = NODES{14,iibox};
                
                
                
                N2 = length(xvalue_rest);
                %Goff = zeros(N1, N2); % off diagonal G, panel 1 as receivers, panel 2 as sources
                Goff = h^2*(1i/4)*besselh(0,k*sqrt((repmat(NODES{06,iibox},1,N2)-xvalue_rest.').^2 ...
                    +(repmat(NODES{07,iibox},1,N2)-yvalue_rest.').^2));
                
                if param.print_flag == 1
                    fprintf(param.fileID,'-------------------\n');
                    fprintf(param.fileID,'Row ID:\n');
                end
                % U*Goff(J,:) ~ Goff
                [J,U,rk] = ID_row(Goff,param.compress_acc);
                %ktot = ktot + rk;
                % store ID for the leaf level
                NODES{08,iibox} = J;
                NODES{09,iibox} = U;
                
                
                if param.print_flag == 1
                    fprintf(param.fileID,'Size of the matrix: %d by %d \n',size(Goff,1),size(Goff,2));
                    errf = norm(U*Goff(J,:) - Goff,'f');
                    nm = norm(Goff,'f');
                    fprintf(param.fileID,'Frob norm is %d\n',errf/nm);
                    fprintf(param.fileID,'Compressed rank is %d\n',rk);
                end
                
                % as the first panel of this level
                % store the compressed GF between son1 and son2 panel
                N1 = length(J); N2 = N1;
                %N2 = length(J2);
                
                skelxp1 = NODES{06,iibox}(J);
                skelyp1 = NODES{07,iibox}(J);
                
                %             skelxp2 = NODES{06,son2}(J2); % skeleton x points from son1
                %             skelyp2 = NODES{07,son2}(J2); % skeleton y points from son2
                
                skelxp2 = NODES{06,iibox+1}(J);
                skelyp2 = NODES{07,iibox+1}(J);
                
                
                Goff12 = h^2*(1i/4)*besselh(0,k*sqrt((repmat(skelxp1,1,N2)-skelxp2.').^2 ...
                    +(repmat(skelyp1,1,N2)-skelyp2.').^2));
                NODES{15,iibox} = Goff12;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % further factorization of compressed G
                [Q,R,rkk] = rank_qr(Goff12,param.compress_acc);
                % store the factorization
                NODES{16,iibox} = Q; NODES{17,iibox} = R;
                
                if param.print_flag == 1
                    fprintf(param.fileID,'-------------------\n');
                    fprintf(param.fileID,'Further factorization of compressed G:\n');
                    fprintf(param.fileID,'Size of the matrix: %d by %d \n',size(Goff12,1),size(Goff12,2));
                    errf = norm(Q*R - Goff12,'f');
                    nm = norm(Goff12,'f');
                    fprintf(param.fileID,'Frob norm is %d\n',errf/nm);
                    fprintf(param.fileID,'Compressed rank is %d\n',rkk);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
            end
            %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
            if param.plot_flag == 1
                plot(NODES{06,iibox} + 1i*NODES{07,iibox},'o');
                plot(NODES{06,iibox}(J) + 1i*NODES{07,iibox}(J),'x');
            end
            %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
        elseif NODES{05,iibox} == 2 % if it has children(non-leaf) node
            
            if level ~= 1 % if non-root level
                
                
                
                if iibox == 2^(level-1) % only compress the first elements
                    
                    for ibox = 2^(level-1):(2^(level)-1)
                        % find the children
                        first = 2^(level-1 + 1); % first on the level below
                        son1 = NODES{4,ibox}(1);
                        son2 = NODES{4,ibox}(2);
                        
                        % store the skeleton points from its children
                        NODES{06,ibox} = [NODES{06,son1}(NODES{08,first});NODES{06,son2}(NODES{08,first})];
                        NODES{07,ibox} = [NODES{07,son1}(NODES{08,first});NODES{07,son2}(NODES{08,first})];
                    end
                    xvalue_rest = NODES{13,iibox};
                    yvalue_rest = NODES{14,iibox};
                    
                    
                    
                    N2 = length(xvalue_rest);
                    Goff = h^2*(1i/4)*besselh(0,k*sqrt((repmat(NODES{06,iibox},1,N2)-xvalue_rest.').^2 ...
                        +(repmat(NODES{07,iibox},1,N2)-yvalue_rest.').^2));
                    
                    if param.print_flag == 1
                        fprintf(param.fileID,'-------------------\n');
                        fprintf(param.fileID,'Row ID:\n');
                    end
                    % U*Goff(J,:) ~ Goff
                    [J,U,rk] = ID_row(Goff,param.compress_acc);
                    %ktot = ktot + rk;
                    % store ID
                    NODES{08,iibox} = J;
                    NODES{09,iibox} = U;
                    
                    if param.print_flag == 1
                        fprintf(param.fileID,'Size of the matrix: %d by %d \n',size(Goff,1),size(Goff,2));
                        errf = norm(U*Goff(J,:) - Goff,'f');
                        nm = norm(Goff,'f');
                        fprintf(param.fileID,'Frob norm is %d\n',errf/nm);
                        fprintf(param.fileID,'Compressed rank is %d\n',rk);
                    end
                    
                    % store the compressed GF between son1 and son2 panel
                    N1 = length(J); N2 = N1;
                    %N2 = length(J2);
                    
                    skelxp1 = NODES{06,iibox}(J);
                    skelyp1 = NODES{07,iibox}(J);
                    
                    %             skelxp2 = NODES{06,son2}(J2); % skeleton x points from son1
                    %             skelyp2 = NODES{07,son2}(J2); % skeleton y points from son2
                    
                    skelxp2 = NODES{06,iibox+1}(J);
                    skelyp2 = NODES{07,iibox+1}(J);
                    
                    
                    Goff12 = h^2*(1i/4)*besselh(0,k*sqrt((repmat(skelxp1,1,N2)-skelxp2.').^2 ...
                        +(repmat(skelyp1,1,N2)-skelyp2.').^2));
                    NODES{15,iibox} = Goff12;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % further factorization of compressed G
                    [Q,R,rkk] = rank_qr(Goff12,param.compress_acc);
                    % store the factorization
                    NODES{16,iibox} = Q; NODES{17,iibox} = R;
                    if param.print_flag == 1
                        fprintf(param.fileID,'-------------------\n');
                        fprintf(param.fileID,'Further factorization of compressed G:\n');
                        fprintf(param.fileID,'Size of the matrix: %d by %d \n',size(Goff12,1),size(Goff12,2));
                        errf = norm(Q*R - Goff12,'f');
                        nm = norm(Goff12,'f');
                        fprintf(param.fileID,'Frob norm is %d\n',errf/nm);
                        fprintf(param.fileID,'Compressed rank is %d\n',rkk);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end
                %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
                if param.plot_flag == 1
                    plot(NODES{06,iibox} + 1i*NODES{07,iibox},'o');
                    plot(NODES{06,iibox}(J) + 1i*NODES{07,iibox}(J),'x');
                end
                %title(strcat(num2str(iibox),'=',num2str(NODES{04,iibox}(1)),'+',num2str(NODES{04,iibox}(2))));
                %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
            end  % compression only for non-root level
            if param.plot_flag == 1
                plot(NODES{06,1} + 1i*NODES{07,1},'o') % root level points
            end
        end
        
    end
    %ktot_last = ktot;
    if param.plot_flag == 1
        hold off;
    end
    %%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%
end
end