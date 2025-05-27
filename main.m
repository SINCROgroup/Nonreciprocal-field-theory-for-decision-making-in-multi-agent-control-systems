%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Agent Based (AB) SIMULATIONS

clc
close
clear

%%%%%%%%%%%%%%%%%%% PARAMETERS OF THE SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% NUMBER OF AGENTS
N = 400;    %Number of Herders
M =400;     %Number of Targets

%%%%%%%%% Parameters of the soft reciprocal repulsion between all the agents
k_rep = 100;    %repulsion strength
sigma = 1;      %particle diameter = repulsion range

D = 1;          % Noise amplitude

%%%%%%% SPACIAL DOMAIN
%%% Lengths are measured in units of sigma
%%% The domain is a squared box of size L centred around the origin
L = 80;

correction = 1;  %%% If correction it is 1, we are considering a periodic domain, and that  we are using the minimum image convention to compute the distances.

%%%%%%% TIME DOMAIN %%%%%%%%
% time is measured in units of \sigma^2/D

dt = .001;          % integration step
time = 300;         % total duration of the simulation
t_settling=200;     % settling time after which we start to save the positions of the agents
frame_spacing = 100; %%% Positions are saved every frame_spacing integration steps

%%%%%%%% Parameters for the repulsion exerted on targets by herders
kt =3;          %repulsion strength
lambda = 2.5;   %repulsion range

%%%%%%%% Parameters for the attraction exerted on herders by targets. If
%%%%%%%% gamma>0 or delta>0 this attraction also incorporates the
%%%%%%%% decion-making of the herders. 

kh = 3;         %attraction strenght
xi = 5;         %attraction range
TD=0;

% if TD (Target Division) is 1, herders only chase a target if in their sensing range there is no other herder closer to that target, enforcing cooperation.
% if TD=0, the herders do not cooperate and can end up chasing the same target

%%%%%%%%%%%%%%%%%%%%%%%% decision-making parameters %%%%%%%%%%%%%%%%%%%%%%
gamma=10;       %parameter tuning selection of the target
delta=2.5;      %parameter tuning the directedness of the chasing position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params=[N,M,k_rep,sigma,D,L,dt,time,t_settling,frame_spacing,kt,lambda,kh,xi,TD];

directoryName="Data_AB";
mkdir(directoryName)  %where the data will be saved


%%  shepherding with circular goal region
% gamma>=0 and delta>=0
AB_radial(params,gamma,delta,directoryName)

%%  shepherding with rectangular goal region (containment or expulsion)
% gamma>=0 and delta>=0 (for containment), or % gamma<=0 and delta<=0 (for expulsion)
AB_rectangular(params,gamma,delta,directoryName)

%%  shepherding with multiple rectangular goal regions (static patterns)
% gamma>=0 and delta>=0 if you want the prescribed goal regions to be
% containment goal regions

goal_regions_pos=[-L/3,0,L/3]; %Centers of the goal regions

AB_static_patterns(params,gamma,delta,goal_regions_pos,directoryName)

%%  shepherding producing travelling patterns
% gamma and delta must have the same sign. 
% if both are positive, you can generate a travelling pattern moving left;
% viceversa it will move right

AB_travelling_patterns(params,gamma,delta,directoryName)

%% Video player AB

x0 = -L/2;
xf = L/2;
y0 = -L/2;
yf = L/2;

%%%% load the data
load("Data_AB/AB_travelling_patterns_g100000_d25000.mat")
t_length=length(H_save(1,1,:));

for it=1:t_length
        figure(1)
        scatter(periodic(squeeze(H_save(:,1,it)),x0,xf),periodic(squeeze(H_save(:,2,it)),x0,xf), MarkerFaceColor='b',SizeData=30)
        hold on
        scatter(periodic(squeeze(T_save(:,1,it)),x0,xf),periodic(squeeze(T_save(:,2,it)),x0,xf), MarkerFaceColor='magenta',SizeData=30)
        hold off
        title(sprintf("%.1f",dt*it*frame_spacing))
        xlim([x0,xf])
        ylim([y0,yf])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Field Equations Simulations

%%%%%%%%%%% DOMAIN
clc
clear
close

L =  80;            % size of the 1d periodic domain centered around the origin
N = 200;            % Number of grid points

dt = .0001;         % Time integration step
T = 3;              % Total integration time
numSteps = T / dt;  %Number of time steps


%%%%%%%%%%%%%%% PARAMETERS
D=5;                % Linear diffusion coefficient
r0=.5;              % Average density 

kt=3;                           % Microscopic coupling for the nonreciprocal repulsion
lambda=2.5;                     % Microscopic interaction range of the nonreciprocal repulsion

kh=3;                           % Microscopic coupling for the nonreciprocal attraction
xi=2.5;                         % Microscopic interaction range of the nonreciprocal attraction

krep=75;                            % Microscopic coupling for the reciprocal repulsion
sigma=1;                            % Microscopic interaction range of the reciprocal repulsion (and unit length)
SRR=(1/3)*krep*(sigma^3)*(r0/D);    % Coarse grained coefficient of the reciprocal repulsion

rot=0.5;                            %initial homogeneous density for targets and herders
roh=0.5;

params=[L;N;dt;T;D;r0;kt;lambda; kh; xi;sigma;SRR;rot;roh];

DirectoryName="DATA_PDE";
mkdir(DirectoryName)  %where the data will be saved

gamma=2.5;
delta=1.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "type" is a string specifying the type of v1 and v2 coupling functions used in the
%simulation of the field equation of the herders' density.

type="main";        % main: v1 and v2 as derived from the microscopic model

% The following 4 choices reproduce the behaviours in Fig 5. v2 is always a
% positive constant

% type="containment"; %containment: v1 is a square wave with  2\pi/L period
% 
% type="expulsion"; %expulsion:  v1 is a square wave with  2\pi/L period,
% but with inverse sign
% 
% type="static_patterns"; %static patterns:  v1 is a square wave with  2\pi/(3L) period
% 
% type="travelling_patterns"; %expulsion:  v1 is a positive constant. This
% specific choice uses slightly different parameters, as discussed in the
% paper. They are specified in "PDE_simulation".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDE_simulation(params,gamma, delta,type, DirectoryName)

%% Video player PDE
clc

%%%% load the data
load("Data_PDE/PDE_main_gamma2500_delta1250.mat")
t_length=length(uH_save(1,:));

for it=1:t_length
    figure(1)
        plot(x,uH_save(it,:),LineWidth=2.2,Color='b')
        hold on
        plot(x,uT_save(it,:),LineWidth=2.2,Color='magenta')
        hold off
        xlim([-(L/2),(L/2)])
        title(sprintf("time %.1f",it*samp_time))
end

%% Fig 4 simulations parameters

%%
%arrays of values for gamma and delta to reproduce the AB part of Fig 4
gamma_array=10.^linspace(-3,1,11);
delta_array=10.^linspace(log10(0.01),log10(0.5),10);


%%
%arrays of values for gamma and delta to reproduce the PDE part of Fig 4 
gamma_array=10.^linspace(-3,log10(2.5),10);
delta_array=2.5*10.^linspace(log10(0.01),log10(0.5),10);






