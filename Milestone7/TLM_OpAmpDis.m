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

g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.05;

% Input parameters for the left side
InputParasL.E0 = 0e7;
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
plotN = 500;

L = 1000e-6*1e2; %u cm
XL = [0,L];
YL = [0,InputParasL.E0];

Nz = 500;
dz = L/(Nz-1);
dt = dz/v_g;
fsync = dt*v_g/dz;

Nt = 200000;
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


beta_spe = .3e-5;
gamma = 1.0;
SPE = 7;

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
    
    beta_r = 0;
    gain_z = gain.*(N - Ntr)./v_g;
    beta_i = (gain_z - alpha)./2;
    beta = beta_r + 1i*beta_i;
    
    InputL(i) = Ef1(t2,InputParasL);
    %InputL(i) = Ef1(t2,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i);
    Er(Nz) = InputR(i);

    exp_det = exp(-1i*dz*beta);

    %Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz - 1);

    %Er(1:Nz-1) = fsync*Er(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz);

    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        Tf = (randn(Nz,1) +1i*randn(Nz,1))*A;
        Tr = (randn(Nz,1) +1i*randn(Nz,1))*A;
    else
        Tf = (ones(Nz,1))*A;
        Tr = (ones(Nz,1))*A;
    end

    EsF = Tf*abs(SPE).*sqrt(N.*1e6);
    EsR = Tr*abs(SPE).*sqrt(N.*1e6);

    Ef = Ef + EsF;
    Er = Er + EsR;


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


    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1)-Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1)-Pr(2:Nz-1));



    if mod(i,plotN)==0
        subplot(4,2,1)
        plot(z*1e6, N, 'b')
        xlabel('z (\mum)')
        %xlim([0 1000])
        ylabel('N')
        ylim([0 4e18])
        title('Carrier Density as a function of z')

        subplot(4,2,2)
        plot(z*1e6, Ef, 'r')
        xlabel('z (\mum)')
        %xlim([0 1000e-6])
        ylabel('Ef')
        ylim([-2e6 2e6])
        title('Stimulated Emission as a function of z')

        subplot(4,1,3)
        plot(time(1:i)*1e12, Nave(1:i), 'g')
        xlabel('Time (ps)')
        ylabel('Nave')
        xlim([0 5000])
        title('Average Carrier Density versus Time')

        subplot(4,1,4);
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12, real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12, real(OutputL),'m');
        xlim([;0,Nt*dt*1e12])
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off

        pause(0.01)
    end

    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;

end


% Compute Fourier transforms of the outputs
fftOutputL = fftshift(fft(InputL));
fftOutputR = fftshift(fft(OutputR));

omega = fftshift(wspace(time));

% Plot the magnitude spectra
figure;
plot(omega/2*pi, 20.*log10(abs(fftOutputL)), 'r', 'LineWidth', 2);
hold on;
plot(omega/2*pi, 20.*log10(abs(fftOutputR)), 'b', 'LineWidth', 2);
hold off;
xlabel('w');
ylabel('Magnitude');
legend('Left Input', 'Right Output');

% Plot the magnitude spectra with logarithmic scale
figure;
plot(omega/2*pi, unwrap(angle(fftOutputL)), 'r', 'LineWidth', 2); % Using semilogy to plot with logarithmic scale
hold on;
plot(omega/2*pi, unwrap(angle(fftOutputR)), 'b', 'LineWidth', 2);
hold off;
xlabel('Frequency (Hz)');
ylabel('Magnitude (log scale)'); % Update ylabel to indicate logarithmic scale
legend('Left Input', 'Right Output');
