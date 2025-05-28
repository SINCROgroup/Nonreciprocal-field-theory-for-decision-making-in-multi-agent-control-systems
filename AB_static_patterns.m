
function AB_static_patterns(params,gamma,delta,goal_regions_pos,DirectoryName)

%%% This function simulates the Agent Based equations in the case the goal
%%% is to contain the targets in three equispaced rectangular regions whose
%%% centres' positions are specified by the array "goal_regions_pos" 


%%%%%%%%%%%%%
% gamma is the parameter tuning the selection of the targets. If gamma=0
% the herders are equally attracted to all the targets. If gamma>0, the
% herders preferentially chases targets with larger distance from the
% closest goal region

% delta is the parameter tuning the trajectory planning during the chasing.
% If delta=0, the herder approach the targets (or the selected target if
% gamma>> 1/xi) from a random direction. If delta >0, the herder apprach
% the target(s) positioning itself at a distance delta at the back of the
% target (at the back wrt the centre of the closest goal region).

%DirectoryName is the directory where the trajectories are saved

%%%%%%% NUMBER OF AGENTS
N = params(1);    %Number of Herders
M =params(2);     %Number of Targets

%%%%%%%%% Parameters of the soft reciprocal repulsion between all the agents
k_rep = params(3);    %repulsion strength
sigma = params(4);      %particle diameter = repulsion range

D = params(5);          % Noise amplitude

%%%%%%% SPACIAL DOMAIN
%%% Lengths are measured in units of sigma
%%% The domain is a squared box of size L centred around the origin
L = params(6);
x0 = -L/2;
xf = L/2;
y0 = -L/2;
yf = L/2;
correction = 1;  %%% If correction it is 1, we are considering a periodic domain, and that  we are using the minimum image convention to compute the distances.

%%%%%%% TIME DOMAIN %%%%%%%%
% time is measured in units of \sigma^2/D

dt = params(7);      % integration step
time = params(8);     % total duration of the simulation
time_steps = round(time/dt);

t_settling=params(9);                         % settling time after which we start to save the positions of the agents
settling_steps=round(t_settling/dt);

frame_spacing = params(10); %%% Positions are saved every frame_spacing integration steps

%%%%%%%% Parameters for the repulsion exerted on targets by herders
kt =params(11);          %repulsion strength
lambda = params(12);   %repulsion range

%%%%%%%% Parameters for the attraction exerted on herders by targets. If
%%%%%%%% gamma>0 or delta>0 this attraction becomes incorporates the
%%%%%%%% decion-making of the herders. Otherwise it is a standard linear
%%%%%%%% attraction

kh = params(13);         %attraction strenght
xi = params(14);         %attraction range
TD = params(15);
% if TD (Target Division) is 1, herders only chase a target if in their sensing range there is no other herder closer to that target, enforcing cooperation.
% if TD=0, the herders do not cooperate and can end up chasing the same target

%%%%%%% DATA SAVING %%%%%%%
save_numbers=round((time_steps-settling_steps)/frame_spacing); % number of saved frames

H_save=zeros(N,2,save_numbers); %arrays where to save herders and targets positions
T_save=zeros(M,2,save_numbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INITIAL POSITIONS generated randomly and uniformly

H0 = initial_pos(N, x0, xf, y0, yf);
T0 = initial_pos(M, x0, xf, y0, yf);


%%% Initialise the arrays storing the deterministic, non-reciprocal
%%% interactions for targets (IT) and herders (IH)
IT = zeros(M,2);
IH = zeros(N,2);

%%% Initialise positions
H=H0;
T=T0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter=1; %%%% counter for saved frames is initialized

for it = 1:time_steps

    % Compute distances

    % herder-target distances
    dHT_x = minimum_image_distance(H(:,1), T(:,1)', L, correction);
    dHT_y = minimum_image_distance(H(:,2), T(:,2)', L, correction);
    dHT = sqrt(dHT_x.^2 + dHT_y.^2);

    % herder-herder distances
    dHH_x = minimum_image_distance(H(:,1), H(:,1)', L, correction);
    dHH_y = minimum_image_distance(H(:,2), H(:,2)', L, correction);
    dHH = sqrt(dHH_x.^2 + dHH_y.^2);
    
    % target-target distances
    dTT_x = minimum_image_distance(T(:,1), T(:,1)', L, correction);
    dTT_y = minimum_image_distance(T(:,2), T(:,2)', L, correction);
    dTT = sqrt(dTT_x.^2 + dTT_y.^2);

    %%% reciprocal soft repulsion

    SRR_pair_HT = repulsion(dHT, -dHT_x, -dHT_y, k_rep, sigma);
    SRR_pair_TT = repulsion(dTT, dTT_x, dTT_y, k_rep, sigma);
    SRR_pair_HH = repulsion(dHH, dHH_x, dHH_y, k_rep, sigma);

    SRR_T = squeeze(sum(SRR_pair_HT, 1)) + squeeze(sum(SRR_pair_TT, 2));
    SRR_H = -squeeze(sum(SRR_pair_HT, 2)) + squeeze(sum(SRR_pair_HH, 2));

    %%% nonreciprocal repulsion exerted on targets by herder
    if kt ~= 0 && lambda ~= 0 %compute repulsion if repulsion strenght and range are nonzero
        IT_pair = repulsion(dHT, -dHT_x, -dHT_y, kt, lambda);
        IT = squeeze(sum(IT_pair, 1));
    end


    %%% nonreciprocal attraction exerted on herders by targets +
    %%% decision-making

    dH0x = minimum_image_distance(H(:,1), goal_regions_pos, L, correction); %minimum image x distance of the herders from the three goal regions
    [~,idxH]=min(abs(dH0x),[],2);                                           % identifies closest goal region to every herder

    for qt=1:N

        dH0x_new(qt,1)=squeeze(dH0x(qt,idxH(qt)));                      %Nx1 array with the x distance of every herder from the closest goal region

    end
    dH0x=dH0x_new;

    %%% same as done above, but for the targets
    dT0x = minimum_image_distance(T(:,1), goal_regions_pos, L, correction);
    [m,idxT]=min(abs(dT0x),[],2);

    for qt=1:M

        dT0x_new(qt,1)=squeeze(dT0x(qt,idxT(qt)));

    end

    dT0x=dT0x_new;

    if kh ~= 0 && xi ~= 0 %compute interaction/decision-making if interaction strenght and range are non-zero

            if TD==1  %enforcing cooperation. Each herder considers the targets within distance xi if there is no other herder closer to the observed targets within its sensing range

                right=zeros(N,M); % each element r_ij: 1 if herder i has to consider target j; 0 viceversa

                for idx=1:N
                    dHT_copy=dHT;
                    mask=zeros(N,M);
                    mask(dHH(:,idx)<=xi,dHT(idx,:)<=xi)=1;
                    dHT_copy(logical(1-mask))=inf;
                    right(idx,:)=(dHT(idx,:)==min(dHT_copy.*mask,[],1)).*(mask(idx,:)==1);
                end
                for i=1:N
                    if sum(right,1)>1
                        sprintf("error")
                        return
                    end
                end

            else  %no cooperation. Each herder considers targets within distance xi
                right=dHT.*(dHT<=xi);
            end
		
            
            right_norm = right .* exp(-gamma .* abs(dH0x));  % the matrix indicating whether a given herder has to consider a given target is now multiplied by the exponential weight in the selection rule for each target

            dx_chasing = minimum_image_distance(H(:,1), T(:,1)' + delta*sign(dT0x)', L, correction);
            dy_chasing = minimum_image_distance(H(:,2), T(:,2)', L, correction);
            IH = attraction(right_norm, dx_chasing, dy_chasing, abs(dT0x), kh, gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %right_label=sum(right,2);
            %right_label=right_label>0;

            %right_norm = right .* exp(-g .* dH0);

            %HT_star=weighted_avg(right_norm, dHT_x, dHT_y, T_rho, g); %H-T^*
            %[cacu,Theta_star]=cartesian_to_polar(minimum_image_distance(H(:,1)-HT_star(:,1),0,L,correction),minimum_image_distance(H(:,2)-HT_star(:,2),0,L,correction));

            %IHx=-a*right_label.*(HT_star(:,1)-delta*lambda*(cos(Theta_star)) );
            %IHy=-a*right_label.*(HT_star(:,2)-delta*lambda*(sin(Theta_star)) );
            %IH=[IHx,IHy];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    % Update positions
    H = H + (IH + SRR_H) * dt + sqrt(2 * D * dt) * randn(N, 2);
    T = T + (IT + SRR_T) * dt + sqrt(2 * D * dt) * randn(M, 2);


    % Save positions
    if (mod(it,frame_spacing)==0)
        figure(1)
        scatter(periodic(H(:,1),x0,xf),periodic(H(:,2),x0,xf), MarkerFaceColor='b',SizeData=30)
        hold on
        scatter(periodic(T(:,1),x0,xf),periodic(T(:,2),x0,xf), MarkerFaceColor='magenta',SizeData=30)
        hold off
        title(sprintf("%.1f",dt*it))
        xlim([x0,xf])
        ylim([y0,yf])

        if (it>settling_steps)

            H_save(:,:,counter)=H;
            T_save(:,:,counter)=T;
            counter=counter+1;
        end
    end

end


filename=sprintf("%s/AB_static_pattern_g%d_d%d",DirectoryName,round(gamma*10000),round(delta*10000));
save(filename,"H_save", "T_save","params","frame_spacing")

end
