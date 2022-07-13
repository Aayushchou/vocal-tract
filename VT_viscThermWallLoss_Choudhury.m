clear all;

% ------------------------------------------------------------------------ %
%                       Vocal Tract: Visco-thermal Losses: VT3
% ------------------------------------------------------------------------ %

% This script consists of the vocal tract model with lip radiation losses, 
% visco-thermal losses and wall damping. The system from VT2 is coupled
% with two forced damped oscillators.

% ------------------------------------------------------------------------- %
%                             Global parameters
% ------------------------------------------------------------------------- %


anim_flag = false;                                                         % Flag for pressure animation: 1 for animation, 0 for no animation
input_type = 'glottal';                                                    % Type of input, either 'glottal' or 'impulse'
vowel = 'e';                                                               % Choice of vowel (a,e,i,o,u)

% ---------- Simulation parameters ----------------- %

anim_speed = 1;                                                            % Speed of animation (larger numbers mean faster)
zoom_order = -2;                                                           % Order of y limits ( from -1*10^order to 1*10^order)
Fs = 44100;                                                                % Sample frequency
Tf = 5;                                                                    % duration

% ----------- Physical parameters ------------------- % 

S0 = 0.00025;                                                              % Area of opening
c = 343;                                                                   % speed of sound

% ---------- Damping parameters -------------------- %

R0 = 100;                                                                  % Resistance
G0 = 100;                                                                  % Conductance
omgp0 = 0;                                                                 % Fundamental of walls 
sigp0 = 8125;                                                              % Damping factor
omgv0 = 0;
sigv0 = 8125; 
rho = 1.225;                                                               % Density of air 
m_wall = 4.76;                                                             % Mass per unit area

% ------------- Glottal waveform parameters --------------- %

depth = 100;                                                               % depth of vibrato
fmod = 0.5;                                                                % Pitch modulation frequency
fmax = 200;                                                                % Maximum frequency
breathiness = 0.0;
% ------------------------------------------------------------------------- %
%                             Derived parameters
% ------------------------------------------------------------------------- %

k = 1/Fs;                                                                  % time step
NF = floor(Tf*Fs);                                                         % number of total samples
L = set_length(vowel);

gamma = c/L;                                                               % gamma parameter
h = gamma*k; N = floor(1/h); h=1/N; lambda = (gamma*k)/h;                  % determine h, N, lambda

if anim_flag == 1                                                          % If animation flag is set, animate
    anim = anim_speed;
else
    anim = NF+1;                                                           % If flag not set, skip animation
end

assert(lambda <= 1);                                                       % Stability condition
assert(R0>=0);                                                             % Assert non-negative resistance
assert(G0>=0);                                                             % Assert non-negative conductance
assert(R0*G0<=1e7);                                                        % Setting an upper bound on damping

% ------------------------------------------------------------------------- %
%                     Cross-sectional area manipulations
% ------------------------------------------------------------------------- %

S = getS_choudhury(vowel);                                                 % Get cross-sectional area function

S = interp1(S(:,1), S(:,2),(0:h:1))';                                      % Interpolating to 1D, normalising between 0 to 1
S(isnan(S)==1) = 1;

S_hat =  S;                                                                % Creating interleaved cross-sectional area approximation
Slhalf = S(1:length(S)-1);

S_hat(2:N) = (Slhalf(2:N) + Slhalf(1:N-1))/2;                              % S and slhalf are related by mu_{x-}

Slhalf(isnan(Slhalf) == 1) = 1;                                            % Setting NaN values to initial cross-sectional area
S_hat(isnan(S_hat) == 1) = 1;                                              % Setting NaN values to initial cross-sectional area

Sav =[S_hat(1); 0.25*(S_hat(3:N+1)+2*S_hat(2:N)+S_hat(1:N-1)); S_hat(N+1)];% Average surface area functions at boundaries and body

% ------------------------------------------------------------------------- %
%                   Loss coefficients for velocity
% ------------------------------------------------------------------------- %

R1 = (2*R0)/((2)+(k*R0));

alpha_v = ((2*rho) - (k*R0))/((2*rho) + (k*R0));                           % alpha loss coefficient for v
beta_v = (2*(k/h))/(2*rho + k*R0);                                         % beta loss coefficient for v

alpha_qv = (2*k*R1)/(2*rho + k*R0);                                        % Alpha q values for velocity

tau_v = ((2)-(k*R0))/((2)+(k*R0));                                         % Loss coefficient tau for determining vhat
xi_v = (k*R0)/((2)+ (k*R0));                                               % Loss coefficient eps for determing vhat

% ------------------------------------------------------------------------- %
%                   Loss coefficients for pressure
% ------------------------------------------------------------------------- %

G1 = (2*G0)/(2 + k*G0);
E = (rho*(gamma^2)*k)/(2 + G0*k);                                          % E parameter

alpha_p = (1 - (E*G0))/(1 + (E*G0));                                       % Alpha loss coefficient for p
beta_p = (((k*rho*(gamma^2)))./((S_hat.*h)*((1+(k*L)))));                  % Beta loss coefficient for pressure 

alpha0_p = (2*E*G0)/(1+(E*G0));                                            % Alpha0 loss coefficient for determining p
alphaq_p = (2*E*G1)/(1+(E*G0));                                            % Alpha q loss coefficient for determining p

xi_p = -(2*G1*k)/(2 + G0*k);                                               % vi loss coefficient for determining p0
eps_p = (2 - G0*k)/(2 + G0*k);                                             % epsilon loss coefficient for determining p0
nyu_p = (G0*k)/(2 + G0*k);                                                 % Nyu loss coefficient for determining p0

% ------------------------------------------------------------------------ %
%                              Lip Radiation
% ------------------------------------------------------------------------ %

lamhalf = 0.5*lambda^2;                                                    % Lambda squared parameter
lam2 = 2*lambda^2;                                                         % Lambda squared parameter
  
alpha0 = (L)/(0.8216*sqrt((S0*S(N))/(pi)));                                % Loss parameter 2 
alpha1 = 1/(2*0.8216^2*gamma);                                             % Loss parameter 1

q1 = (alpha0*lambda^2*h*((3*S_hat(N) - S_hat(N-1))/(2*S_hat(N))));         % Boundary coefficient                                                     
q2 = (alpha1*lambda^2*h*((3*S_hat(N) - S_hat(N-1))/(2*k*S_hat(N))));       % Boundary coefficient

r1 = lam2 / (1+q1+q2);                                                     % defines radiation 
r2 = -(1+q1-q2)/(1+q1+q2);                                                 % defines radiation

% ------------------------------------------------------------------------ %
%                   Wall Loss damping for pressure
% ------------------------------------------------------------------------ %

eps = c * sqrt((2*rho)/(m_wall))*((pi/S0)^(0.25));                         % Coupling coefficient 

% ------------ oscillator coefficients ------------- %

WLp0 = (1+(sigp0*k));                                                      % Coefficient for n+1 terms
WLp1 = (2 - ((k^2)*(omgp0^2)));                                            % Coefficient for n terms
WLp2 = (1-(sigp0*(k)));                                                    % Coefficient for n-1 terms
WLp3 = (0.5*eps*S_hat.^(0.25))*k./(Sav);                                   % Coupling parameter
WLp4 = (eps*S_hat.^(0.25))*k./(2);                                         % Coupling parameter


wcoefp0 = WLp1./WLp0;                                                      % Calculate coefficient
wcoefp1 = WLp2./WLp0;                                                      % Calculate coefficient
wcoefp2 = gamma*WLp4./WLp0;                                                % Calculate coefficient
wcoefp3 = (1./(1+WLp3.*WLp4));                                             % Coefficient for whole scheme 
% ------------ Coupling coefficients -------------- %


pcoup0 = (WLp3*WLp1)/WLp0;                                                 % Coupling coefficient
pcoup1 = (WLp3*2)/WLp0;                                                    % Coupling coefficient
pcoup2 = ((WLp3.*WLp4)/(WLp0));                                            % Coupling coefficient

% ------------------------------------------------------------------------ %
%                   Wall Loss damping for velocity
% ------------------------------------------------------------------------ %

eps = c * sqrt((2*rho)/(m_wall))*((pi/S0)^(0.25));                         % Coupling coefficient 

% ------------ oscillator coefficients ------------- %

WLv0 = (1+(sigv0*k));                                                      % Coefficient for n+1 terms
WLv1 = (2 - ((k^2)*(omgv0^2)));                                            % Coefficient for n terms
WLv2 = (1-(sigv0*(k)));                                                    % Coefficient for n-1 terms
WLv3 = (0.5*eps*Slhalf.^(0.25))*k./(Slhalf);                               % Coupling parameter
WLv4 = (eps*Slhalf.^(0.25))*k./(2);                                        % Coupling parameter


wcoefv0 = WLv1./WLv0;                                                      % Calculate coefficient
wcoefv1 = WLv2./WLv0;                                                      % Calculate coefficient
wcoefv2 = gamma*WLv4./WLv0;                                                % Calculate coefficient
wcoefv3 = (1./(1+WLv3.*WLv4));                                               % Coefficient for whole scheme
% ------------ Coupling coefficients -------------- %


vcoup0 = (WLv3*WLv1)/WLv0;                                                 % Coupling coefficient                                               
vcoup1 = (WLv3*2)/WLv0;                                                    % Coupling coefficient
vcoup2 = ((WLv3.*WLv4)/(WLv0));                                            % Coupling coefficient

% ------------------------------------------------------------------------- %
%                             Defining input
% ------------------------------------------------------------------------- %

if strcmp(input_type,'impulse')    
    x = zeros(NF, 1);                                                      % Impulse input
    x(1) = 1;
    
elseif strcmp(input_type,'glottal')
    t = 0:1/Fs:Tf;                                                         % Time vector
    f_in = (depth/2)*cos(2*pi*fmod*t') + fmax-(depth/2);                   % frequency vector
    phase_in = cumsum(f_in/Fs);                                            % integral of frequency
    x = sin(2*pi*phase_in);                                                % generating pitch-varying sinusoid
    epsilon = 0.8;                                                             
    x = abs(x-epsilon)/(1-epsilon);                                        % Applying abs function to get glottal waveform shape 
    x(1:end-1) = x(1:end-1).*(1-breathiness) +...                          % incorporating breathiness
        sin(x(1:end-1).*breathiness.*diff(rand(length(x),1)));
else
    
    error('Please select valid input type: glottal, impulse')
end


% ------------------------------------------------------------------------- %
%                           System matrices
% ------------------------------------------------------------------------- %

d1 = ones(N,N+1);
d1(:,1:2:N+1) = -1;
D1 = spdiags(d1,0:1,N,N+1);                                                % 22x21 matrix mapping pressure to velocity
D2 = -D1';                                                                 % 21x22 matrix mapping velocity to pressure 
D2(N+1,N) = 0;
% ------------------------------------------------------------------------- %
%                             Initialising
% ------------------------------------------------------------------------- %

p = zeros(N+1, 1); p1 = zeros(N+1,1); p2 = zeros(N+1,1);                   % Initialising vectors for pressure

p0 = zeros(N+1,1); p01 = zeros(N+1,1); ptil = zeros(N+1,1); ...            % Initialising phat,p0           
ptil1 = zeros(N+1,1);

v = zeros(N,1); v1 = zeros(N,1); vpr = zeros(N,1); vpr1 = zeros(N,1);      % Initialising velocity vectors
v2 = zeros(N,1);
wp = zeros(N+1,1); wp1 = zeros(N+1,1); wp2 = zeros(N+1,1);                 % initialising w vectors
wv = zeros(N,1); wv1 = zeros(N,1); wv2 = zeros(N,1);
out = zeros(NF,1);                                                         % Initialing output vector 

% ------------------------------------------------------------------------- %
%                             Update scheme
% ------------------------------------------------------------------------- %
tic;
for n= 1:NF  
            
    v = wcoefv3.*(alpha_v.*v - beta_v.*(D1*p)...                                     % Updating velocity
        + alpha_qv.*vpr ...
        - vcoup0.*wv1 + vcoup1.*wv2 + vcoup2.*v2);   
    
    vpr = tau_v*vpr + xi_v*(v + v1);                                       % Updating v'
    
    v(1) = x(n);                                                           % Exciting input 
            
    p = wcoefp3.*(alpha_p*p - ...                                                    % Updating pressure
        beta_p.*(D2*(v.*Slhalf))...
        + alpha0_p*p0 + alphaq_p*ptil ...
        - pcoup0.*wp1 + pcoup1.*wp2 + pcoup2.*p2); 
    
    p0 = eps_p*p01 + nyu_p*(p+p1)...                                       % Update p0
        + xi_p*ptil1;
    
    ptil = -ptil1 + ...
        (p + p1 - p0 - p01);                                               % Update ptilde
    
    p(N) = r1*p1(N-1) + r2*p2(N);                                          % Applying radiational losses
        
    wp = wcoefp0*wp1 - wcoefp1*wp2 + wcoefp2.*(p-p2);                      % Calculating wall vibrations 
    
    wv = wcoefv0*wv1 - wcoefv1*wv2 + wcoefv2.*(v-v2);                      % Calculating wall vibrations 
    
    if mod(n,anim) == 0                                                    % Drawing the pressure outputs  
        if n == 1        
            figure
            plt1 = plot(p);
            xlim([0 23])                                                   % Limiting axis
            ylim([-5*(10^(zoom_order)) 5*(10^(zoom_order))])
        else
            set(plt1, 'ydata', p)
            drawnow
        end       
    end
    
    out(n) = p(N);                                                         % Reading output pressure    
    v2 = v1; v1 = v; vpr1 = vpr; p01 = p0; ptil1 = ptil;p2=p1; p1 = p;
    
end
toc 
% ------------------------------------------------------------------------- %
%                             Plotting outputs
% ------------------------------------------------------------------------- %

f_ax = (0:NF-1)*(Fs/NF);                                                   % Frequency axis
semilogy(f_ax,abs(fft(out)))                                               % Plotting impulse response
xlim([0 5000])                                                             % Limitting x-axis
ylim([10e-1 10e4])
hold on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title...
('Formant structure comparison /e/: Vocal tract visco-thermal losses')

out = out/max(abs(out));                                                   % Normalise output
soundsc(out(200:end),Fs);                                                  % Output sound


function [L] = set_length(vowel)
        if strcmp(vowel,'a')
            L = 0.1746;
        elseif strcmp(vowel,'e')
            L = 0.1667;
        elseif strcmp(vowel,'i')
            L = 0.1667;
        elseif strcmp(vowel,'o')
            L = 0.1746;
        elseif strcmp(vowel,'u')
            L = 0.1746;
        end
end


