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
InputParasL.E0 = 1e5;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 5e-13;
InputParasL.phi = 0;

% Input parameters for the right side 
InputParasR = 0;

% Group velocity and wavelength
n_g = 3.5;
vg = c_c/n_g*1e2; %TWM cm/s group velocity
Lambda = 1550e-9;

% Simulation parameters
plotN = 10;

L = 1000e-6*1e2; %cm
XL = [0,L];
YL = [0,InputParasL.E0];

Nz =500;
dz = L/(Nz-1);
dt = dz/vg;
fsync = dt*vg/dz;

Nt = floor(2*Nz);
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

Ef1 = @SourceFct;
ErN = @SourceFct;

t = 0;
time(1) = t;

InputL(1) = Ef1(t,InputParasL);
InputR(1) = ErN(t,0);

OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

Ef(1) = InputL(1);
Er(Nz) = InputR(1);

figure('name', 'Fields')
subplot(3,1,1)
plot(z*10000,real(Er),'b');
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z\mum')
ylabel('E_r')
hold off
subplot(3,1,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputL),'r--'); 
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputR),'b--'); 
xlabel('time(ps)')
ylabel('E')

hold off

for i = 2:Nt
    t = dt*(i-1);
    time(i) = t;
    
    InputL(i) = Ef1(t,InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i);
    Er(Nz) = InputR(i);

    Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*Er(2:Nz);

    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

    if mod(i,plotN)==0
        subplot(3,1,1)
        plot(z*10000,real(Ef), 'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,1,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')

        hold off
        subplot(3,1,3);
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12, real(OutputR),'g');
        plot(time*1e12,real(InputR),'b');
        plot(time*1e12, real(OutputL),'m');
        xlim([;0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off
        pause(0.01)
    end



end
