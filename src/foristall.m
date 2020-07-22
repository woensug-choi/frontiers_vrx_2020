
% Forristall
A = 42.0;
B = 63.0;

f = logspace(-6,3,1000);
%f = linspace(0,1,20000);
z = 10;
VV = [1, 2, 4, 8, 16];

figure(1);
clf();
figure(2);
clf();

lstr = {};
q = zeros(length(VV),1);
wind_vars = zeros(length(VV),1);
for ii = 1:length(VV)
    vz = VV(ii);
    % variance as a function of mean speed
    w = 1;
    sigma = 3*w*(0.00076 * vz^2 +0.0304*vz);

    % Normalized frequency
    fstar = f*z/vz;
    % Non-dimensional spectrum
    ss = A*fstar./((1+B*fstar).^(5/3));
    % Dimensional spectrum
    s = sigma^2*A*fstar./(f.*(1+B*fstar).^(5/3));
    %s = sigma^2*(A*z/v)./((1+B*z/v*f).^(5/3));
    % Used to check to see if the area under the curve is indeed equal to
    % sigma^2
    q(ii) = trapz(f,s);
    wind_vars(ii) = sigma^2;
    %qq = trapz(ss,fstar)
    
    figure(1)
    subplot(211)  % Non-dimensional
    %semilogx(f,ss)
    loglog(f,ss)
    hold on
    subplot(212) % Dimensional
    %semilogx(f,s)
    loglog(f,s)
    %plot(f,s)
    hold on
    lstr{ii} = sprintf('$v_{10}$=%4.1f m/s',vz);%,q(ii));
    
    % Find bandwidth, cuttoff freq
    %dc = s(1);
    %II = find(s<dc/2);
    %fc1 = f(II(1));
    %xline(fc1,'r')
    %fc2 = (2^(3/5)-1)/B*v/z;
    %xline(fc2,'g')
    %subplot(211)
    %xline(fc2)
    %yline(dc/2)
    
    
    figure(2)
    semilogx(fstar,ss);
    hold on
    

    
end

figure(1)
ax1 = subplot(212);
ylabel('Dimensional Spectra: $S(f)$ [(m/s)$^2$/Hz]','interpreter','latex')
grid on
xlabel('Dimensional Frequency, $f$ [Hz]','interpreter','latex')
ax2 = subplot(211);
grid on
%ylabel('Normalized Spectra: $f \, S(f)/\sigma^2$ [n/a]','interpreter','latex')
ylabel('Normalized Spectra: $\tilde{S}(f)$ [n/a]','interpreter','latex')
legend(lstr,'interpreter','latex','location','south')
linkaxes([ax1, ax2],'x')
xlim(10.^[-5,1])

figure(2)
grid on
%ylabel('Normalized Spectra: $f \, S(f)/\sigma^2$ [n/a]','interpreter','latex')  
ylabel('Normalized Spectra: $\tilde{S}(f^*)$ [n/a]','interpreter','latex')
xlabel('Normalized Frequency, $f^*$ [n/a]','interpreter','latex')
xlim(10.^[-5,3])

% Checking that I calculated the peak freq. correctly by looking at where the
% gradient goes to xer.
% figure(3);
% semilogx(fstar,gradient(ss))
% xlabel('f^*')
% ylabel('dS/df*')
% title('Checking that I calculated the peak freq. correctly')
% grid on


%% Generate timeseries


% Frequency limit -hig enough that S(omega) ~ 0
omega_u = 1.0*2*pi;

% Number of freq samples
N = 5000;
% Find N for T0 = 60
%T0 = 1; %60;
%N = round(T0 * 2*omega_u / (2*pi));
% Freq. sample size
domega = omega_u/N*2;

% Time series
% Final time - one full period of generated series
T0 = 1 * 2*pi/domega;

% Number of time steps
M = 10000;
% time step
dt = T0/M;
% Time vector
tt1 = (0:dt:(M-1)*dt)';
% Random seed
s = rng(42);
% Generate phases for all runs
phases = rand(N,1)*2*pi;

figure(4);
clf()
figure(5);
clf()
figure(6);
clf()

% Values for mean velocity
VZ = [1, 2, 4, 8, 16];
VZ = [1, 4, 16];
lstr = {};
hh = [];
for ii = 1:length(VZ)
    % Forristall parameters
    % Height
    z = 10;
    % Mean velocity at z
    vz = VZ(ii); %10;

    % Initial time series 
    yy1 = zeros(length(tt1),1);
    % Loop through, starting with 2 so that An for n = 0 is zero
    for n = 2:N-1
        % Sample the PSD
        omegan = n*domega;
        An = sqrt(1/(2*pi)*forristall_spectra(omegan/(2*pi),z,vz)*domega);
        yy1 = yy1 + sqrt(2)*An *cos(omegan*tt1+phases(n+1));
    end
    
    % Normalized timeseries
    % Theory
    w = 1;
    sigma = 3*w*(0.00076 * vz^2 +0.0304*vz);

    lstr{ii} = sprintf('$v_{10}$=%4.1f m/s, $\\sigma$=%4.2f m/s',vz,sigma);%,q(ii));
    sigma_est = std(yy1)
    
    figure(4)
%     if ii == 1
%         yyaxis left
%     else
%         yyaxis right
%     end
    subplot(5,1,2:3)
    plot(tt1,yy1/sigma);
    hold on
    
    % Non-normalized
    %figure(6)
    subplot(5,1,4:5)
    hh(ii) = plot(tt1,yy1);
    hold on
    
    % Sanity check
    fprintf('Mean vel.: %f\n',vz)
    fprintf('Mean: %f\n',mean(yy1))
    fprintf('Var: %f\n',var(yy1))
    
    % Calculate PSD to compare with generating function.
    Nx = length(yy1);
    nsc = floor(Nx/4.5);
    nov = floor(nsc/2);
    nff = max(256,2^nextpow2(nsc));
    fs1 = 1/median(diff(tt1));
    [Sig1,f1]=pwelch(yy1,hamming(nsc),nov,nff,fs1);
    
    % Dimensional spectra
    figure(5)
    loglog(f1,Sig1)
    hold on
    
    % Normalized spectra
end
legend(lstr,'interpreter','latex')
grid on


figure(4)
subplot(5,1,2:3)
TT = 2*60;
xlim([0,TT])
grid on;
ylabel('$v_g(t)/\sigma$ [n/a]','interpreter','latex')
%figure(6)
subplot(5,1,4:5)
xlim([0,TT])
grid on;
ylabel('$v_g(t)$ [m/s]','interpreter','latex')
xlabel('Time [s]')
hLegend = subplot(5,1,1)
posLegend = get(hLegend,'Position');
%leg = legend(hlegend,hh,lstr,'interpreter','latex');
leg = legend(hh,lstr,'interpreter','latex');
axis(hLegend,'off');
set(leg,'Position',posLegend);

%figure(7)
%plot(f1,'.')