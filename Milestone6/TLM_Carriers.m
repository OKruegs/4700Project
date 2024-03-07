% Setting default properties for plots
set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'Defaultaxeslinewidth', 2)
set(0, 'DefaultFigureWindowStyle', 'docked')

% Constants and parameters
c_c = 299792458;              % Speed of light
c_eps_0 = 8.8542149e-12;      % Vacuum permittivity (F/m)
c_eps_0_cm = c_eps_0 / 100;   % Vacuum permittivity in F/cm
c_mu_0 = 1 / c_eps_0 / c_c^2; % Magnetic constant
c_q = 1.60217653e-19;         % Elementary charge
c_hb = 1.05457266913e-34;     % Reduced Planck's constant
c_h = c_hb * 2 * pi;          % Planck's constant


g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.1;

beta_r = 0; %rotation
beta_i = 0; %gain

% Input parameters for the left side
InputParasL.E0 = 5e7;       %Ef Scale
InputParasL.we =  0;%   20e12; %Frequency 
InputParasL.t0 = 5e-12;     %Delay/Time
InputParasL.wg = 5e-13;%10e-6;     %Width 
InputParasL.phi = 0;        %Phase 

% Input parameters for the right side 
InputParasR.E0 = 0;%4e7;       %Er Scale
InputParasR.we = 0;            %Frequency
InputParasR.t0 = 0;%2e-12;     %Delay/Time
InputParasR.wg = 0;%4e-13;     %Width 
InputParasR.phi = 0;           %Phase  

% Group velocity and wavelength
Ntr = 1e18;
n_g = 3.5;
v_g = c_c/n_g*1e2; % TWM cm/s group velocity
Lambda = 1550e-9; %cm
f0 = c_c/Lambda;

% Simulation parameters
plotN = 20;           %Speed
L = 1000e-6 * 2e2;    % Length of the system in cm
XL = [0, L];          % X-axis limits
YL = [-1e7, 1e7]; % Y-axis limits
Nz = 1000;             % Number of spatial points
dz = L / (Nz - 1);    % Spatial step size delta x
dt = dz / v_g;         % Temporal step size delta t
fsync = dt * v_g / dz; % Synchronization factor

Nt = floor(2 * Nz);   % Total number of time steps
tmax = Nt * dt;
t_L = dt * Nz;        % Time to travel length

% Initialization
z = linspace(0, L, Nz).'; % Spatial grid
time = nan(1, Nt);        % Time array
InputL = nan(1, Nt);      % Left input wave
InputR = nan(1, Nt);      % Right input wave
OutputL = nan(1, Nt);     % Left output wave
OutputR = nan(1, Nt);     % Right output wave

%Milestone 6 parameters
N = ones(size(z))*Ntr;
Nave(1) = mean(N);
% if GenGifs
%     system(['rm' gitFile]);
% end
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

Ef = zeros(size(z)); % Forward wave
Er = zeros(size(z)); % Reverse wave

Pf = zeros(size(z));
Pr = zeros(size(z));

Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

Ef1 = @SourceFct2; % Function for left input source
ErN = @SourceFct2; % Function for right input source 

t = 0;
time(1) = t;

InputL(1) = Ef1(t, InputParasL);
InputR(1) = ErN(t, InputParasR);

OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

Ef(1) = InputL(1);
Er(Nz) = InputR(1);

% Define parameters for the grating
kappa0 = 0;       % Maximum value of kappa
kappaStart = 1/3;   % Starting position of the grating (as a fraction of the total length)
kappaStop = 2/3;    % Ending position of the grating (as a fraction of the total length)

% Initialize kappa with maximum value everywhere
kappa = kappa0 * ones(size(z));

% Set kappa to 0 outside the grating region
kappa(z < L * kappaStart) = 0;
kappa(z > L * kappaStop) = 0;

% Plot initial fields
figure('name', 'Fields')
subplot(3, 1, 1)
plot(z * 10000, real(Er), 'b');
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3, 1, 2)
plot(z * 10000, real(Er), 'b');
xlabel('z\mum')
ylabel('E_r')
hold off
subplot(3, 1, 3)
plot(time * 1e12, real(InputL), 'r'); hold on
plot(time * 1e12, real(OutputL), 'r--');
plot(time * 1e12, real(InputR), 'b'); hold on
plot(time * 1e12, real(OutputR), 'b--');
xlabel('time(ps)')
ylabel('E')
hold off

% Main simulation loop
for i = 2:Nt
    t = dt * (i - 1);
    time(i) = t;
    

    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, InputParasR);

     % Update boundary reflections
     RL = 0.9i; % Reflection coefficient for the left boundary
     RR = 0.9i; % Reflection coefficient for the right boundary
     Ef(1) = InputL(i);% + RL * Er(1);   %Reflection
     Er(Nz) = InputR(i);% + RR * Ef(Nz); %Reflection

%     beta = ones(size(z))*(beta_r+1i*beta_i);
%     exp_det = exp(-1i*dz*beta);

%     Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz - 1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);   %Update equations
%     Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(1:Nz-1).*Ef(1:Nz-1);   


%     Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz - 1);
%     Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz);

    % Update reflected outputs
    OutputR(i) = Ef(Nz);% * (1 - RR); %Reflection
    OutputL(i) = Er(1);% * (1 - RL);  %Reflection

%     Pf(1) = 0;
%     Pf(Nz) = 0;
%     Pr(1) = 0;
%     Pr(Nz) = 0;
%     Cw0 = -LGamma + 1i*Lw0;
% 
%     Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
%     Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
%     Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
%     Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);
% 
%     Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
%     Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));

S = (abs(Ef).^2 +abs(Er).^2).*EtoP*1e-6;

if t < Ion || t > Ioff
    I_injv = I_off;
else
    I_injv = I_on;
end

Stim = gain.*(N-Ntr).*S;
N = (N + dt*(I_injv/eVol - Stim))./(1+ dt/taun);
Nave(i) = mean(N);

    % Plotting every 'plotN' steps
    if mod(i, plotN) == 0
%         subplot(3, 1, 1)
%         plot(z * 10000, real(Ef), 'r'); hold on
%         plot(z * 10000, imag(Ef), 'r--'); hold off
%         xlim(XL * 1e4)
%         ylim(YL)
%         xlabel('z(\mum)')
%         ylabel('E_f')
%         legend('\Re', '\Im')
%         hold off
% 
%         subplot(3, 1, 2)
%         plot(z * 10000, real(Er), 'b'); hold on
%         plot(z * 10000, imag(Er), 'b--'); hold off
%         xlim(XL * 1e4)
%         ylim(YL)
%         xlabel('z(\mum)')
%         ylabel('E_r')
%         legend('\Re', '\Im')
%         hold off
% 
%         subplot(3, 1, 3);
%         plot(time * 1e12, real(InputL), 'r'); hold on
%         plot(time * 1e12, real(OutputR), 'g');
%         plot(time * 1e12, real(InputR), 'b');
%         plot(time * 1e12, real(OutputL), 'm');
%         xlim([0, Nt * dt * 1e12])
%         ylim(YL)
%         xlabel('time(ps)')
%         ylabel('E')
%         legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
%         hold off
% 
%         pause(0.01)


        % Plot N as a function of z
  
        subplot(3,1,1);
        plot(z * 10000, N(:, end)); % Plot the final timestep
        xlabel('z (\mum)');
        ylabel('N');
        title('N as a function of z');

        % Plot S as a function of z
        subplot(3,1,2);
        plot(z * 10000, S(:, end)); % Plot the final timestep
        xlabel('z (\mum)');
        ylabel('S');
        title('S as a function of z');

        % Plot Nave as a function of time
        subplot(3,1,3);
        plot(time * 1e12, Nave);
        xlabel('Time (ps)');
        ylabel('N_{ave}');
        title('N_{ave} as a function of time');

        pause(0.01);
    end



    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;

end



