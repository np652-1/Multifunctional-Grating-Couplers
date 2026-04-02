function [Lvnew,Nlfv,Nv,Npv,Lhnew,Nlfh,Nh,Nph] = identical_grid_info(Lv,Lh,h)
% generate uniform grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to generate identical panels on the leaf level given a given domain [Lv,Lh] with spacing h,
% this function automatically extends the computational domain [Lv,Lh] to [Lvnew,Lhnew].
% It does so in a way to keep [Lvnew,Lhnew] as close to [Lv,Lh] as possible
% while making sure the panel on the leaf level has about 100 points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs:
%     Lv: user input vertical length
%     Lh: user input horizontal length
%     h:  user input grid spacing

% This function is written in a way that the new domain will never be
% smaller than the input domain and h is unchanged, for the convenience of
% the user


%Outputs:
%     Lvnew % new vertical domain length
%     Nlfv  % leaf panel vertical point number
%     Nv    % total vertical points number
%     Npv   % total vertical panel number
% 
%     Lhnew % new horizontal domain length
%     Nlfh  % leaf panel horizontal point number
%     Nh    % total horizontal points number
%     Nph   % total horizontal panel number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    powmax = 29; % max number of divisions. Change to larger number if the computation is unusually large
    Lvnew_opt = inf; % a big number for initializaiton
    Lhnew_opt = inf; % a big number for initializaiton


    for Nlfv = 8:15
    %for Nlfv = 5:9
        for Nlfh = 8:15
        %for Nlfh = 10:18

            Npv = 1; % vertical panel number
            Nph = 1; % horizontal panel number

        for i = 1:powmax

            ifmorev = (Nlfv*Npv - 1)*h < Lv;
            ifmoreh = (Nlfh*Nph - 1)*h < Lh;

            if ifmorev
                Npv = Npv*2;
            end

            if ifmoreh
                Nph = Nph*2;
            end

            if ifmorev == 0 && ifmoreh == 0
                break;
            end
        end

    if i == 29
        error('More than 2^%d = %d panel on one side. Too many panels.', powmax,2^powmax);
    end

            Nv = Npv*Nlfv;
            Nh = Nph*Nlfh;

            %fprintf('-----------\n'); fprintf('-----------\n');

            Lvnew = (Nv - 1)*h;
            Lhnew = (Nh - 1)*h;



            if Lhnew < Lhnew_opt            
                Lhnew_opt = Lhnew;
                Nh_opt = Nh;
                Nlfh_opt = Nlfh;
                Nph_opt = Nph;
            end


            if Lvnew < Lvnew_opt
                Lvnew_opt = Lvnew;
                Nv_opt = Nv;
                Nlfv_opt = Nlfv;
                Npv_opt = Npv;
            end
            %fprintf('-----------\n'); fprintf('-----------\n');
        end
    end


    Lvnew = Lvnew_opt; % new vertical domain length
    Nlfv = Nlfv_opt;   % leaf panel vertical point number
    Nv = Nv_opt;       % total vertical points number
    Npv = Npv_opt;     % total vertical panel number

    Lhnew = Lhnew_opt;% new horizontal domain length
    Nlfh = Nlfh_opt;  % leaf panel horizontal point number
    Nh = Nh_opt;       % total horizontal points number
    Nph = Nph_opt;     % total horizontal panel number
end

