function NODES = inverse_build(NODES,k,h,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Follow Algorithm 1 in the Gopal Martinsson paper.
%
%       This function takes an upward pass from leaf level to construct the
%       inverse of the integral operator (I + B*G). 
%       It is achieved by computing the scattering matrix S by first
%       inverting each diagonal block of (I + B*G) on each leaf panel.
%       Then merge the S matrices of two children to form S of their
%       parent.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%       Inputs:
%       NODES structure with compressed Green's function stored (local_compression)
%                            spatial Chi distribution BB stored (get_geometry)
%       k: wave number
%       h: grid spacing
%       param (direct solver parameters):
%             param.quad_order: quadrature order of the direct solver
%             (Only support 2 or 4 for the current version)
%             param.D0: correction for 4th order quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%       Outputs:
%       NODES structure with inverse and S matrices computed and stored
%                    NODES{11,ibox}   X
%                    NODES{12,ibox}   S
%       They are stored to hierarchically in order to apply the inverse in O(NlogN) time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Upward pass to construct all scattering matrices
    nlevels = NODES{01,end};
    for level = nlevels:-1:1
        first = 2^(level-1); % first of the level
        U = NODES{09,first};
        
        if level < nlevels
            first = 2^(level-1 + 1); % first of the level below
            Goff12 = NODES{15,first};
            Q = NODES{16,first}; R = NODES{17,first};
            rk = size(Q,2);
            %norm(Q*R - Goff12)  
        end
        
        for ibox =   2^(level-1):(2^level-1) % leaf level, solve the problem on each panel
            
            if NODES{05,ibox} == 0 % if node node
                
                %first = 2^(level-1); % first of the level
                
                x1  = NODES{06,ibox}; N1 = length(x1);
                y1  = NODES{07,ibox};
                %J1 =  NODES{08,ibox};
                %U =  NODES{09,first};
                B1 =  NODES{10,ibox};
                
                if ibox == 2^(level-1) % only do it once
                    Gdiag = h^2*(1i/4)*besselh(0,k*sqrt((repmat(x1,1,N1)-x1.').^2 ...
                                                       +(repmat(y1,1,N1)-y1.').^2));
                    switch param.quad_order
                        case 2
                            Gdiag(isnan(Gdiag))= 0;
                        case 4
                            Gdiag(isnan(Gdiag))= (1i/4)*param.D0*h^2;
                        otherwise
                            error('Invalid quadrature order number.');
                    end
                end

                % scattering matrix
                X = inv(eye(N1) + B1*Gdiag); 
                NODES{11,ibox} = X;
                S = (U.')*X*(B1*U); 
                NODES{12,ibox} = S;
                %T = eye(N1) + B1.*Gdiag;
                %S = (U1.')*(T\(B1.*U1));
                %fprintf(fileID,'%d\n',S);

        elseif NODES{5,ibox} == 2 %

                %first = 2^(level-1 + 1); % first of the level below
                %level = nlevels-1;

                %for ibox =   2^(level-1):(2^level-1) % parent node
                %end


                son1 = NODES{4,ibox}(1);
                son2 = NODES{4,ibox}(2);



                N1 = size(Goff12,1); N2 = N1;

                S1 = NODES{12,son1}; % scattering matrix from son1
                S2 = NODES{12,son2}; % scattering matrix from son2



%                  XX = inv([   eye(N1)   ,  S1*Goff12         ;... 
%                                                                  ...
%                          S2*(Goff12.')     ,  eye(N2)   ]);
%             

%                 X = eye(N1+N2) - ...
%                               [ S1*Q,zeros(N1,rk);zeros(N2,rk),S2*R.']...
%                               *inv([eye(rk), R*S2*R.' ; Q.'*S1*Q, eye(rk)])...
%                               *[zeros(rk,N1),R; Q.',zeros(rk,N2)];
                X = eye(N1+N2) - ...
                              [ S1*Q,zeros(N1,rk);zeros(N2,rk),S2*R.']...
                              *([eye(rk), R*S2*R.' ; Q.'*S1*Q, eye(rk)]\[zeros(rk,N1),R; Q.',zeros(rk,N2)]);
                %norm(XX - X)  
                NODES{11,ibox} = X;

                if level ~= 1 % if not the root level
                    %first = 2^(level-1); % first of the level
                    %U = NODES{09,first};
                    % combine scattering matrix from children 
                    
                    Scom = [    S1      , zeros(N1,N2)      ;    ... 
                                                                 ...
                           zeros(N2,N1) ,      S2          ];
                    % scattering matrix
                    S = (U.')*X*Scom*U;
                    NODES{12,ibox} = S;
                    %fprintf(fileID,'%d\n',S);
                end

            end
        end

    end
end
