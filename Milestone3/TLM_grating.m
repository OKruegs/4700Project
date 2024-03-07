%Comments:What each block is doing

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

%Real part of beta causes rotating field
%Imaginary part of beta causes simple gain

beta_r = 0;
beta_i = 0; %nothing seems to change

%beta_r = 0;
%beta_i = 8; %gain occurs as expected (Ran with simple and modulated)

%beta_r = 80;
%beta_i = 0; %no gain but modulation and stability falls

%beta_r = 80;
%beta_i = 8;  %Increased gain and modulation on simple gaussian

% Input parameters for the left side
InputParasL.E0 = 300;       %Ef Scale
InputParasL.we = 0;%0; %Frequency 
InputParasL.t0 = 3e-12;     %Delay/Time
InputParasL.wg = 5e-13;     %Width 
InputParasL.phi = 0;        %Phase 

% Input parameters for the right side 
InputParasR.E0 = 0;       %Er Scale
InputParasR.we = 0;            %Frequency
InputParasR.t0 = 2e-12;     %Delay/Time
InputParasR.wg = 5e-13;     %Width 
InputParasR.phi = 0;           %Phase  

% Group velocity and wavelength
n_g = 3.5;
vg = c_c / n_g * 1e2; % Group velocity (cm/s)
Lambda = 1550e-9;     % Wavelength

% Simulation parameters
plotN = 20;           %Speed
L = 1000e-6 * 1e2;    % Length of the system in cm
XL = [0, L];          % X-axis limits
YL = [-InputParasL.E0, InputParasL.E0]; % Y-axis limits
Nz = 300;             % Number of spatial points
dz = L / (Nz - 1);    % Spatial step size delta x
dt = dz / vg;         % Temporal step size delta t
fsync = dt * vg / dz; % Synchronization factor

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

Ef = zeros(size(z)); % Forward wave
Er = zeros(size(z)); % Reverse wave

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
kappa0 = 300;       % Maximum value of kappa
kappaStart = 1/3;   % Starting position of the grating (as a fraction of the total length)
kappaStop = 2/3;    % Ending position of the grating (as a fraction of the total length)

% Initialize kappa with maximum value everywhere
kappa = kappa0 * ones(size(z));

% Set kappa to 0 outside the grating region
kappa(z < L * kappaStart) = 0;
kappa(z > L * kappaStop) = 0;

Eprevf = Ef; % Store initial values of Ef

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
%     RL = 0.9i; % Reflection coefficient for the left boundary
%     RR = 0.9i; % Reflection coefficient for the right boundary
    Ef(1) = InputL(i);% + RL * Er(1);   %Reflection
    Er(Nz) = InputR(i);% + RR * Ef(Nz); %Reflection

    beta = ones(size(z))*(beta_r+1i*beta_i);
    exp_det = exp(-1i*dz*beta);

    Ef_temp = Ef;

    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz - 1) + 1i*dz*kappa(2:Nz).*Er(2:Nz);   %Update equations
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i*dz*kappa(1:Nz-1).*Ef_temp(1:Nz-1);   

    Ef_temp = Ef;

    % Update reflected outputs
    OutputR(i) = Ef(Nz);% * (1 - RR); %Reflection
    OutputL(i) = Er(1);% * (1 - RL);  %Reflection

    % Plotting every 'plotN' steps
    if mod(i, plotN) == 0
        subplot(4, 1, 1)
        plot(z * 10000, real(Ef), 'r'); hold on
        plot(z * 10000, imag(Ef), 'r--'); hold off
        xlim(XL * 1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re', '\Im')
        hold off
        subplot(4, 1, 2)
        plot(z * 10000, real(Er), 'b'); hold on
        plot(z * 10000, imag(Er), 'b--'); hold off
        xlim(XL * 1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re', '\Im')

        hold off
        subplot(4, 1, 3);
        plot(time * 1e12, real(InputL), 'r'); hold on
        plot(time * 1e12, real(OutputR), 'g');
        plot(time * 1e12, real(InputR), 'b');
        plot(time * 1e12, real(OutputL), 'm');
        xlim([0, Nt * dt * 1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('E')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off

        subplot(4, 1, 4);
        phase = angle(Ef);
        plot(z * 10000, phase, 'k');
        xlim(XL * 1e4)
        xlabel('z (\mum)')
        ylabel('Phase')
        hold off
        

        pause(0.01)
    end


end

% Plot kappa as a function of z
figure;
plot(z*10000, kappa); % Multiply by 10000 to convert from meters to micrometers for better visualization
xlabel('z (\mum)');
ylabel('\kappa');
title('Grating Profile (\kappa as a function of z)');
grid on;

% Compute Fourier transforms of the outputs
fftOutputL = fft(InputL);
fftOutputR = fft(OutputR);

% Shift zero frequency component to the center
fftOutputL_shifted = fftshift(fftOutputL);
fftOutputR_shifted = fftshift(fftOutputR);

% Define the frequency range
N = length(OutputL);
fs = 1 / (time(2) - time(1)); % Sampling frequency
frequencies = linspace(-fs/2, fs/2, N); % Frequency range

% Plot the magnitude spectra
figure;
plot(frequencies, abs(fftOutputL_shifted), 'r', 'LineWidth', 2);
hold on;
plot(frequencies, abs(fftOutputR_shifted), 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Left and Right Outputs');
legend('Left Input', 'Right Output');
grid on;

% Plot the magnitude spectra with logarithmic scale
figure;
semilogy(frequencies, abs(fftOutputL_shifted), 'r', 'LineWidth', 2); % Using semilogy to plot with logarithmic scale
hold on;
semilogy(frequencies, abs(fftOutputR_shifted), 'b', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude (log scale)'); % Update ylabel to indicate logarithmic scale
title('Frequency Spectrum of Left and Right Outputs (Logarithmic Scale)');
legend('Left Input', 'Right Output');
grid on;


