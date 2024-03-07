% Setting default properties for plots
set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth',2);
set(0, 'Defaultaxeslinewidth',2)

set(0,'DefaultFigureWindowStyle','docked')

% Constants and parameters
c_c = 299792458;             % Speed of light
c_eps_0 = 8.8542149e-12;     % Vacuum permittivity (F/m)
c_eps_0_cm = c_eps_0 / 100;   % Vacuum permittivity (F/cm)
c_mu_0 = 1 / c_eps_0 / c_c^2; % Magnetic constant
c_q = 1.60217653e-19;         % Elementary charge
c_hb = 1.05457266913e-34;     % Reduced Planck's constant hbar
c_h = c_hb * 2 * pi;          % Planck's constant

% Input parameters for the left side
InputParasL.E0 = 1e7;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13;
InputParasL.phi = 0;

Ntr = 1e18;
n_g = 3.5;
v_g = c_c/n_g*1e2; % TWM cm/s group velocity
Lambda = 1550e-9; % cm
f0 = c_c/Lambda;

% Simulation parameters
plotN = 50;

L = 1000e-6*1e2; %u cm
XL = [0,L];
YL = [0,InputParasL.E0];

Nz = 500;
dz = L/(Nz-1);
dt = dz/v_g;
fsync = dt*v_g/dz;

Nt = 100000;
%Nt = floor(2*Nz);
tmax = Nt*dt;
t_L = dt*Nz; %time to travel length

z = linspace(0,L,Nz).'; %Nz points, Nz-1 segments
time = nan(1,Nt);
InputL = nan(1,Nt);
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

Ef = zeros(size(z));
Er = zeros(size(z));

Ef1 = @SourceFct3;
ErN = @SourceFct3;

t = 0;
time(1) = t;

InputL(1) = Ef1(t,InputParasL);
InputR(1) = ErN(t,0);

OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

Ef(1) = InputL(1);
Er(Nz) = InputR(1);

%Milestone 6 parameters
N = ones(size(z))*Ntr;
Nave(1) = mean(N);

gain = v_g*2.5e-16; 
eVol = 1.5e-10*c_q;
Ion = 0.25e-9;
Ioff = 3e-9;
I_off = 0.024;
I_on = 0.1;
taun = 1e-9;
Zg = sqrt(c_mu_0/c_eps_0)/n_g;
EtoP = 1/(Zg*f0*v_g*1e-2*c_hb);
alpha = 0;

p = 0;
for i = 2:Nt

    t = dt*(i-1);
    time(i) = t;

    if mod(i,10000) == 0
        p = 0;
    else
        p = p+1;
        t2 = dt*(p-1);
    end
    
    InputL(i) = Ef1(t2,InputParasL);
    %InputL(i) = Ef1(t2,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i);
    Er(Nz) = InputR(i);

    Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*Er(2:Nz);

    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

    % Carrier code
    S = (abs(Ef).^2 +abs(Er).^2).*EtoP*1e-6;

    Stim = gain.*(N-Ntr).*S;
    if t < Ion || t > Ioff
        I_injv = I_off;
    else
        I_injv = I_on;
    end
    N = (N + dt*(I_injv/eVol - Stim))./(1+ dt/taun);
    Nave(i) = mean(N);

    if mod(i,plotN)==0
        subplot(3,1,1)
        plot(z*1e6, N, 'b')
        xlabel('z (\mum)')
        xlim([0 1000])
        ylabel('N')
%         ylim([0.5e18 5e18])
        title('Carrier Density as a function of z')

        subplot(3,1,2)
        plot(z*1e6, S, 'r')
        xlabel('z (\mum)')
        %xlim([0 1000e-6])
        ylabel('S')
        title('Stimulated Emission as a function of z')

        subplot(3,1,3)
        plot(time(1:i)*1e12, Nave(1:i), 'g')
        xlabel('Time (ps)')
        ylabel('Nave')
        xlim([0 5000])
        title('Average Carrier Density versus Time')

        pause(0.01)
    end

end