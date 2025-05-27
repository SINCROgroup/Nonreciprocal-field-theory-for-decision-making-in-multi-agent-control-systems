function PDE_simulation(params,gamma,delta,behav_type,DirectoryName)

% integrates in time two coupled density fields in a 1d periodic boundary
% domain.

%The partial differential equations are integrated using a pseudospectral code combined with an operator splitting
%technique, which allows to accurately treat the linear part of the spatial operator (the linear diffusion). The
%non-linear parts of the spatial operator (non reciprocal interactions and reciprocal SR repulsion) are treated as source
%terms: at every time-step, they are evaluated in real space using the values from the previous step. Then they are
%transformed into the Fourier space, where they are considered as source terms and integrated using a fourth-order
%Runge-Kutta time integration scheme


%%%%%%%%%%% DOMAIN

L =  params(1);          % size of the 1d periodic domain centered around the origin
N = params(2);          % Number of grid points

x = linspace(-(L/2), (L/2), N+1);           % Points in the space grid
x = x(1:N);                             
k = [0:round(N/2)-1, -round(N/2):-1] * (2*pi/L);  % wave vector

dt =params(3);       % Time integration step
T = params(4);          % Total integration time

%%%%%%%%%%%%%%% PARAMETERS
D=params(5);             % Linear diffusion coefficient
r0=params(6);            % Average density

kt=params(7);                           % Microscopic coupling for the nonreciprocal repulsion
lambda=params(8);                       % Microscopic interaction range of the nonreciprocal repulsion
R=kt*(1/3)*lambda^3*(r0/D);             % Coarse grained coefficient of the nonreciprocal repulsion

kh=params(9);                           % Microscopic coupling for the nonreciprocal attraction
xi=params(10);      % Microscopic interaction range of the nonreciprocal attraction

sigma=params(11);                            % Microscopic interaction range of the reciprocal repulsion (and unit length)
SRR=params(12);    % Coarse grained coefficient of the soft reciprocal repulsion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now select the type of v1 and v2 functions. Also specifies the
%%%%%% sampling time at which the densities will be saved

switch behav_type

    case "main"
        samp_time=.25;
        frame_rate=round(samp_time/dt);

        v1_transient=@(x) ( delta*(1-2*gamma*x).*(x-xi)+delta*(x+xi)+gamma*x.*(xi^2-x.^2)+(2/3)*gamma*(x.^3) );
        v1=((2*delta*xi+(2/3)*gamma*(xi^3))*((abs(x)>=xi).*(abs(x)<=(L/2)-xi))+v1_transient(abs(x)).*(abs(x)<xi)+v1_transient((L/2)-abs(x)).*(abs(x)>(L/2)-xi)).*sign(x)*(sigma/D)*r0;

        v1=v1*kh;

        v2_transient=@(x) ( delta*(1-gamma*x).*(xi^2-(x.^2)) +(2/3)*(delta*gamma+1)*(xi^3)-(2/3)*(gamma*x).*(xi^3-(x.^3))+(gamma/2)*(xi^4-(x.^4)) );
        v2=( (2/3) *(delta*gamma+1)*(xi^3).*((abs(x)>=xi).*(abs(x)<(L/2)-xi))+ (abs(x)<xi).*v2_transient(abs(x))+(abs(x)>=(L/2)-xi).*v2_transient((L/2)-abs(x))) *(r0/D);

        v2=v2*kh;

    case "containment"
        samp_time=.25;
        frame_rate=round(samp_time/dt);

        v1=(2*delta*xi+(2/3)*gamma*(xi^3))*(sigma/D)*r0*square(pi*x/L);
        v2=(2/3) *(delta*gamma+1)*(xi^3)*(r0/D);
        v1=v1*kh;
        v2=v2*kh;

    case "expulsion"
        samp_time=.25;
        frame_rate=round(samp_time/dt);

        v1=-(2*delta*xi+(2/3)*gamma*(xi^3))*(sigma/D)*r0*square(pi*x/L);
        v2=(2/3) *(delta*gamma+1)*(xi^3)*(r0/D);
        v1=v1*kh;
        v2=v2*kh;

    case "static_patterns"
        samp_time=.25;
        frame_rate=round(samp_time/dt);

        v1=(2*delta*xi+(2/3)*gamma*(xi^3))*(sigma/D)*r0*square(6*pi*x/L);
        v2=(2/3) *(delta*gamma+1)*(xi^3)*(r0/D);
        v1=v1*kh;
        v2=v2*kh;

    case "travelling_patterns"
        kt=20;
        R=kt*(1/3)*lambda^3*(r0/D);

        dt = .001;       % Passo temporale
        T = 750;          % Tempo finale della simulazione
        samp_time=2.5;
        frame_rate=round(samp_time/dt);

        v1=(2*delta*xi+(2/3)*gamma*(xi^3))*(sigma/D)*r0*ones(size(x));
        v2=(2/3) *(delta*gamma+1)*(xi^3)*(r0/D);
        v1=v1*kh;
        v2=v2*kh;

end

numSteps = T / dt; %Number of time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% INITIAL CONDITIONS

rot=params(13);
roh=params(14);

pert1=.01*randn(size(x));
pert2=.01*randn(size(x));

uH=roh*ones(1,N)+pert1;
uT=rot*ones(1,N)+pert2;

%%%% SAVED ARRAYS FOR THE DENSITIES
uH_save=zeros(T/samp_time +1,N);
uT_save=zeros(T/samp_time +1,N);
uH_save(1,:)=uH;
uT_save(1,:)=uT;
counter=2;

%%%%%%%%%%%%%%% OPERATORS AND Runge-Kutta 4 TIME STEPS

grad=(1i*k);
dt_rk4=[0,dt/2,dt/2,dt];


%%%%%%%%%%%%% TIME INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for step = 1:numSteps

    uH_hat=fft(uH);
    uT_hat=fft(uT);
    qH=zeros(4,N); %% per RK4
    qT=zeros(4,N); %% per RK4

    uH_hat_RK=fft(uH);
    uT_hat_RK=fft(uT);


    for i=1:4

        if i>1

            uH_hat_RK=exp(-(k.^2) *dt_rk4(i)).*(uH_hat_RK+dt_rk4(i)*qH(i-1,:));
            uT_hat_RK=exp(-(k.^2) *dt_rk4(i)).*(uT_hat_RK+dt_rk4(i)*qT(i-1,:));

        end


        grad_uT_hat=grad.*uT_hat_RK;
        grad_uH_hat=grad.*uH_hat_RK;

        grad_uH=ifft(grad_uH_hat,'symmetric');
        grad_uT=ifft(grad_uT_hat,'symmetric');

        uT_RK=ifft(uT_hat_RK,'symmetric');
        uH_RK=ifft(uH_hat_RK,'symmetric');

        uT_grad_uT_hat=fft(uT_RK.*grad_uT);
        uT_grad_uH_hat=fft(uT_RK.*grad_uH);
        v2_uH_grad_uT_hat=fft(v2.*uH_RK.*grad_uT);
        uH_grad_uT_hat=fft(uH_RK.*grad_uT);
        uH_grad_uH_hat=fft(uH_RK.*grad_uH);
        v1_uHuT_hat=fft(v1.*uH_RK.*uT_RK);

        qT(i,:)= exp((k.^2) *dt_rk4(i)).*grad.*(SRR*uT_grad_uT_hat+((SRR+R)*uT_grad_uH_hat));
        qH(i,:)= exp((k.^2) *dt_rk4(i)).*grad.*(SRR*uH_grad_uH_hat+((SRR)*uH_grad_uT_hat -v2_uH_grad_uT_hat)-v1_uHuT_hat);


    end


    uH_hat=exp(-(k.^2) *dt).*(uH_hat+(dt/6)*[1,2,2,1]*qH );
    uT_hat=exp(-(k.^2) *dt).*(uT_hat+(dt/6)*[1,2,2,1]*qT );


    uT=real(ifft(uT_hat,'symmetric'));
    uH=real(ifft(uH_hat,'symmetric'));


    if mod(step,frame_rate)==0

        figure(1)
        plot(x,uH,LineWidth=2.2,Color='b')
        hold on
        plot(x,uT,LineWidth=2.2,Color='magenta')
        hold off
        xlim([-(L/2),(L/2)])
        title(sprintf("time %.1f",step*dt))

        uH_save(counter,:)=uH;
        uT_save(counter,:)=uT;
        counter=counter+1;

    end

end

%%%%%%%%%%%%%%%% SAVE
filename=sprintf("%s/PDE_%s_gamma%d_delta%d",DirectoryName,behav_type,round(gamma*1000),round(delta*1000));
save(filename,"uT_save", "uH_save","x","params","samp_time")

end



