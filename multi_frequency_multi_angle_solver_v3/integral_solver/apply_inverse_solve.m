function [qq,varargout] = apply_inverse_solve(NODES,ff,Npv,Nph)
% Apply inverse as a preconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Follow Algorithm 2 in the Gopal Martinsson paper.
        
%       This function applies the formed inverse to the given (source) vector ff
%       in O(N*log N) time.

%       It does so by first doing an upward pass from leaf to root to
%       compute all necessary incoming expansion by combining results from
%       children to parents.
%       Then it does a downward pass to compute all outgoing expansion by 
%       inherit results from parents to children. 

%       Eventualy the solution qq is assembled in the leaf level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Inputs:
%             NODES structure with stored local inverses X and scattering
%             matrices S.                                                 (inverse_build)
%                            
%             k: wave number
%             h: grid spacing
%             ff: vectorized source vector f = -Einc*B, ff = f(:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Outputs:
%       NODES structure with incoming, outoing expansion, final solution stored
%                    NODES{20,ibox}   store source vector ff on each leaf nodes
%                    NODES{21,ibox}   r = X*f only for leaf nodes
%                    NODES{22,ibox}   rt = (U.')*r
%                    NODES{23,ibox}   wt, only for nonroot nodes
%                    NODES{24,ibox}   qt, only for nonroot nodes



%                    NODES{25,ibox}   q, SOLUTION on each panel!
%                                     Tt is stored only if we choose two outputs 
%                                     [qq,NODES] = apply_inverse_solve(NODES,k,h,ff)
%
%                                     For qq = apply_inverse_solve(NODES,k,h,ff)
%                                     only the assembled solution qq = inv(I + BG)*ff is returned.
%                                     This return format is convenient for
%                                     applying the inverse as a
%                                     preconditioner.
    

    Nh = NODES{2,1}(2);
    Nv = NODES{2,1}(4);
    % loop over leaf nodes, get source
    for ibox = 1:size(NODES,2)
        if NODES{05,ibox} == 0 % if leaf node
            index = NODES{26,ibox};
            NODES{20,ibox} = ff(index);
        end
    end

    nlevels = NODES{01,end};
    %%%% apply_inverse starts
    % upward pass
    for level = nlevels:-1:1
        
        first = 2^(level-1); % first of the level
        U = NODES{09,first};

        if level < nlevels
            first = 2^(level-1 + 1); % first of the level below
            Goff12 = NODES{15,first};
        end
        
        if level == 1
            first = 2^(level-1 + 1); % first of the level below
            Q = NODES{16,first}; R = NODES{17,first};
            %rk = size(Q,2);
        end
        
        for ibox = 2^(level-1):(2^level-1)
            %ibox = 2^(level-1);
            %leaf node
            if NODES{05,ibox} == 0 % if leaf node
                
                X = NODES{11,ibox}; f = NODES{20,ibox};
                r = X*f; 
                NODES{21,ibox} = r;
                
                rt = (U.')*r;
                NODES{22,ibox} = rt;
                %fprintf(fileID,'%d\n',r);
                %fprintf(fileID,'%d\n',rt);
            end


            if (NODES{05,ibox} == 2) && level ~= 1 % if nonleaf node and nonroot
                
                %first = 2^(level-1); % first of the level 
                
                son1 = NODES{4,ibox}(1);
                son2 = NODES{4,ibox}(2); 

                %U = NODES{09,first};
                X = NODES{11,ibox};
                rt1 = NODES{22,son1};
                rt2 = NODES{22,son2};
                rt = (U.')*X*[rt1;rt2];
                NODES{22,ibox} = rt;   
                %fprintf(fileID,'%d\n',r);
                %fprintf(fileID,'%d\n',rt);
            elseif level == 1 % root level
                
                
                son1 = NODES{4,ibox}(1);
                son2 = NODES{4,ibox}(2); 
                X = NODES{11,ibox}; 
                rt1 = NODES{22,son1};
                rt2 = NODES{22,son2};

                
                N1 = size(Goff12,1); N2 = N1;
                
%                 wt = [ zeros(N1,N1), Goff12; ...
%                      (Goff12.')   , zeros(N2,N2)]...
%                      * X * [rt1;rt2];  % top-level solve
%                 wt1 = wt(1:N1); wt2 = wt((N1+1):(N1+N2));

                 temp = X * [rt1;rt2];
                 wt1 = Q*R*temp((N1+1):(N1+N2));
                 wt2 = R.'*Q.'*temp(1:N1);
                 NODES{23,son1} = wt1;
                 NODES{23,son2} = wt2;
                
            end        
        end
    end

    % downward pass up to level above the leaf level
    for level = 2:(nlevels-1)
            first = 2^(level-1); % first of the level 
            U = NODES{09,first};
            if level < nlevels
                first = 2^(level-1 + 1); % first of the level below
                Goff12 = NODES{15,first};
                
                Q = NODES{16,first}; R = NODES{17,first};
            %rk = size(Q,2);
            end
            
        for ibox = 2^(level-1):(2^level-1)
            
            X = NODES{11,ibox}; 
            wt = NODES{23,ibox}; 

            son1 = NODES{4,ibox}(1);
            son2 = NODES{4,ibox}(2); 
            rt1 = NODES{22,son1};
            rt2 = NODES{22,son2};
            S1 = NODES{12,son1}; 
            S2 = NODES{12,son2};
            
            temp =  U*wt;
%             
%             qt = X*([rt1;rt2]- [S1, zeros(size(S1,1),size(S2,2));...
%                                 zeros(size(S2,1),size(S1,2)), S2]*U*wt);
            
            qt = X*([rt1;rt2] - [S1*temp(1:size(S1,2));S2*temp(size(S1,2)+1:end)]);

            N1 = size(Goff12,1); N2 = N1;
                                           
%             wt = [ zeros(N1,N1), Goff12; ...
%                  (Goff12.')   , zeros(N2,N2)]*qt + U*wt;

            wt = [Q*R*qt((N1+1):(N1+N2));R.'*Q.'*qt(1:N1)] + temp;

            NODES{23,son1} = wt(1:N1);
            NODES{23,son2} = wt((N1+1):(N1+N2));

        end
    end


    % build solution q
    %qq = zeros(Nv*Nh,1);
    %     figure;
    %     hold on;
    level = nlevels; % leaf level: construct solution
    first = 2^(level-1); % first of the level 
    U = NODES{09,first};
    for ibox = 2^(level-1):(2^level-1)
        
        
        r = NODES{21,ibox};
        X = NODES{11,ibox};
        B = NODES{10,ibox};
        wt = NODES{23,ibox};
        %q = r - X*(B.*U)*wt; % old, B.*U is expensive
        
        %n = length(B);
        %B = spdiags(B,0,n,n);
        q = r - X*B*U*wt;
        
        
        %q = reshape(q,[12,10]);
        q = reshape(q,[Nv/Npv,Nh/Nph]);
        NODES{25,ibox} = q;
%         index = NODES{26,ibox};
%         qqq(index) = q;
    end
    
    xcutnum = log2(Nph);
    
    for level = (nlevels-1):-1:1
        for ibox = 2^(level-1):(2^level-1)
            son1 = NODES{4,ibox}(1);
            son2 = NODES{4,ibox}(2); 
            
            if level >= (xcutnum + 1)
                NODES{25,ibox} = [ NODES{25,son2}; NODES{25,son1}];
            end
            
            if level < (xcutnum + 1)
                NODES{25,ibox} = [ NODES{25,son1}, NODES{25,son2}];
            end
%             NODES{25,son1} = [];  NODES{25,son2} = [];
%         figure;
%         imagesc(abs(NODES{25,ibox}));
        end
    end
    
	
    %assemble q
    qq = NODES{25,1}(:); % root

%     figure;
%     imagesc(reshape(abs(qq),Nv,Nh));
% 
%     figure;
%     imagesc(reshape(abs(qq-qqq),Nv,Nh));
% 
%     if ~isequal(qq,qqq)
%         error('Assemble is wrong');
%     end
    
    if nargout == 2
        varargout{1} = NODES;
    end
end

