clear all;

% ------------------------------------------------------------------------ %
%                       Vocal Tract: Lip Radiation Losses: VT0
% ------------------------------------------------------------------------ %

% This script consists of the vocal tract model with lip radiation losses.
% It solves websters equation for a variable cross-sectional area. Special
% boundary conditions are applied in order to represent radiation at the
% lip boundary. The glottal waveform is represented using a clipped sine
% wave.

% ------------------------------------------------------------------------ %
%                               Global Parameters
% ------------------------------------------------------------------------ %

% --------------------- Options ------------------- %

input_type = 'glottal';                                                    % Type of input: 'glottal' or 'impulse'
vowel = 'e';                                                               % Choice of vowel (a,e,o,u,i)

% --------------- Physical and simulation parameters ------ %

Fs = 44100;                                                                % Sample rate
L = 0.17;                                                                  % Length of vocal tract
c = 343;                                                                   % Speed of sound
Tf = 5;                                                                    % Duration of simulation
s0 = 0.00025;                                                              % Area of vocal tract opening

% ------------- Glottal waveform parameters --------------- %

depth = 100;                                                               % depth of vibrato
fmod = 0.5;                                                                % Pitch modulation frequency
fmax = 200;                                                                % Maximum frequency
breathiness = 0.0;
% ------------------------------------------------------------------------ %
%                               Derived Parameters
% ------------------------------------------------------------------------ %

S = getS_choudhury(vowel);                                                 % Get cross-sectional area function 

k = 1/Fs;                                                                  % Sample period
Nf = floor(Tf*Fs);                                                         % Duration of simulation in samples

gamma = c/L;                                                               % gamma parameter
h = gamma*k; N = floor(1/h); h=1/N; lambda = gamma*k/h;                    % determine h, N, lambda

S = interp1(S(:,1), S(:,2), (0:h:1))';                                     % interpolate area function pairs
S(isnan(S)==1) = 1;                                                        % Set initial values to 1

Sav = [S(1); 0.25*(S(3:N+1)+2*S(2:N)+S(1:N-1)); S(N+1)];                   % Average surface area functions at boundaries and body
Sr = 1.5*S(N+1) - 0.5*S(N);

assert(lambda <= 1);                                                       % Stability condition

% ------------------------------------------------------------------------ %
%                            Input Waveform
% ------------------------------------------------------------------------ %

if strcmp(input_type,'impulse')
    x = zeros(Nf, 1);                                                      % Impulse input
    x(1) = 1;
    
elseif strcmp(input_type,'glottal')
    t = 0:1/Fs:Tf;                                                         % Time vector
    f_in = (depth/2)*cos(2*pi*fmod*t') + fmax-(depth/2);                   % frequency vector
    phase_in = cumsum(f_in/Fs);                                            % integral of frequency
    x = sin(2*pi*phase_in);                                                % generating pitch-varying sinusoid
    epsilon = 0.8;
    x = abs(x-epsilon)/1-epsilon;                                          % Applying abs function to get glottal waveform shape
    x(1:end-1) = x(1:end-1).*(1-breathiness) +...                          % Incorporating effects of breathiness
        sin(x(1:end-1).*breathiness.*diff(rand(length(x),1)));
end

% ------------------------------------------------------------------------ %
%                              Lip Radiation
% ------------------------------------------------------------------------ %

lamhalf = 0.5*lambda^2;                                                    % Lambda squared parameter
lam2 = 2*lambda^2;                                                         % Lambda squared parameter

alpha0 = (L)/(0.8216*sqrt(s0*S(N+1)/pi));                                  % Loss parameter 2
alpha1 = 1/(2*0.8216^2*gamma);                                             % Loss parameter 1

q1 = (alpha0*lambda^2*h*((3*S(N+1) - S(N))/(2*Sav(N+1))));                 % Boundary coefficient
q2 = (alpha1*lambda^2*h*((3*S(N+1) - S(N))/(2*k*Sav(N+1))));               % Boundary coefficient

r1 = lam2 / (1+q1+q2);                                                     % defines radiation
r2 = -(1+q1-q2)/(1+q1+q2);                                                 % defines radiation

g1 = -(k^2*gamma^2/(h/S(1))) * (3*S(1)-S(2));                              % gain

% ------------------------------------------------------------------------ %
%                        Psi Update Coefficients
% ------------------------------------------------------------------------ %

coef2 = lamhalf*((S(2:N)+S(3:N+1))./Sav(2:N));                             % coefficient for n,l-1 term
coef1 = lamhalf*((S(2:N)+S(1:N-1))./Sav(2:N));                             % coefficient for n,l+1 term
coef0 = 2*(1-lambda^2);                                                    % coefficient fir n, l term

% ------------------------------------------------------------------------ %
%                      Finite Difference Scheme
% ------------------------------------------------------------------------ %

psi = zeros(N+1, 1); psi1 = zeros(N+1, 1); psi2 = zeros(N+1, 1);           % Initialising psi vectors
out = zeros(Nf, 1);                                                        % Initialise output vector

tic;                                                                       % Checking efficiency
for n=1:Nf
    
    psi(2:N) = coef0*psi1(2:N)+coef1.*psi1(1:N-1) - psi2(2:N)...           % Main finite difference update scheme
        + coef2.*psi1(3:N+1);
    
    psi(1) = coef0*psi1(1) + lam2*psi1(2) - psi2(1)+ g1*x(n);              % Boundary condition left
    
    psi(N+1) = r1*psi1(N) + r2*psi2(N+1);                                  % Boundary condition right
    
    out(n) = Fs*(psi(N+1)-psi1(N+1));                                      % Generating output
    
    psi2 = psi1; psi1 = psi;                                               % Updating variables
    
end
toc

% ------------------------------------------------------------------------ %
%                              Plotting Results
% ------------------------------------------------------------------------ %

f_ax = (0:Nf-1)*(Fs/Nf);%*(Fs/length(1:Fs/2+1));                           % Frequency axis
semilogy(f_ax,abs(fft(out)))                                               % Plotting impulse response
xlim([0 10000])                                                            % Limiting x-axis
ylim([10e-1 10e4])                                                         % Limiting y-axis
hold on;
xlabel('Frequency (Hz)');
ylabel('Pressure (dB)');
title...
    ('Formant structure comparison /e/: Vocal tract with and without wall damping')

out = out/max(abs(out));                                                   % Normalise
soundsc(out,Fs);                                                           % Output sound


