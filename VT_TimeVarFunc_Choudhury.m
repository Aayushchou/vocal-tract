

% ------------------------------------------------------------------------ %
%                       Vocal Tract: Time-Varying function
% ------------------------------------------------------------------------ %

% This script has the same functionality as VT4, although it is a function
% that is called by the choir_Choudhury function. It has the added
% functionality of turning MIDI notes to frequency.

% ------------------------------------------------------------------------- %
%                             Global parameters
% ------------------------------------------------------------------------- %

function [out] = VT_TimeVarFunc_Choudhury(opts,sim,inp)



% ---------- Options and damping ----------------- %

anim_flag = opts.anim;                                                     % Flag for pressure animation: 1 for animation, 0 for no animation
input_type = opts.type;
Fs = opts.SR;                                                              % Sample frequency
Tf = opts.Tf;                                                              % duration
anim_speed = opts.animSpeed;                                                % Speed of animation (larger numbers mean faster)
zoom_order = opts.zoom;                                                     % Order of y limits ( from -1*10^order to 1*10^order)
R0 = opts.R0;                                                              % Resistance
G0 = opts.G0;                                                              % Conductance
vowel = opts.vowel;                                                        % Choice of vowel (a,e,i,o,u)


% ----------- Physical parameters ------------------- % 

S0 = 0.00025;                                                              % Area of opening
c = 343;                                                                   % speed of sound
L = sim.L;                                                                 % Length of vocal tract
rho = 1.225;                                                               % Air density


% ------------- Glottal waveform parameters --------------- %

MIDI = inp.MIDI;                                                           % Midi input 
depth = inp.depth;                                                         % depth of vibrato
fmod = inp.fmod;                                                           % Pitch modulation frequency

f0 = 2^((MIDI-69)/12)*440;                                                 % Fundamental frequency
fmax = f0+depth;                                                           % Get maximum frequency
breathiness = inp.breath;                                                  % Breathiness control

% ------------- Vowel change parameters --------------- %

squish = 2;                                                                % How gradually one vowel changes to another
shift = 0.4;                                                               % The position in time of the switch

% ------------------------------------------------------------------------- %
%                             Derived parameters
% ------------------------------------------------------------------------- %

k = 1/Fs;                                                                  % time step
NF = floor(Tf*Fs);                                                         % number of total samples
L = set_length('o');



gamma = c/L;                                                               % gamma parameter
h = gamma*k; N = floor(1/h); h=1/N; lambda = (gamma*k)/h;                  % determine h, N, lambda

interp = 0.5*tanh((squish*(1:NF).*k)-(shift*NF*k)) + 0.5;                   % Interpolation function
interp = interp';


if anim_flag == 1                                                          % If animation flag is set, animate
    anim = anim_speed;
else
    anim = NF+1;                                                           % If flag not set, skip animation
end

assert(lambda <= 1);                                                       % Stability condition

% ------------------------------------------------------------------------- %
%                     Cross-sectional area manipulations
% ------------------------------------------------------------------------- %

[S1, S2] = getAreas(vowel);                                                % Get cross-sectional areas

S1 = interp1(S1(:,1), S1(:,2),(0:h:1))';                                   % Interpolating to 1D, normalising between 0 to 1
S1(isnan(S1)==1) = 1;

S2 = interp1(S2(:,1), S2(:,2),(0:h:1))';                                   % Interpolating to 1D, normalising between 0 to 1
S2(isnan(S2)==1) = 1;

S_hat1 =  S1(1:length(S1));                                                % Creating interleaved cross-sectional area approximation
Slhalf1 = S1(1:length(S1)-1);

S_hat1(2:N) = (Slhalf1(2:N) + Slhalf1(1:N-1))/2;

Slhalf1(isnan(Slhalf1) == 1) = 1;                                          % Setting NaN values to initial cross-sectional area
S_hat1(isnan(S_hat1) == 1) = 1;                                            % Setting NaN values to initial cross-sectional area

S_hat2 =  S2(1:length(S2));                                                % Creating interleaved cross-sectional area approximation
Slhalf2 = S2(1:length(S2)-1);

S_hat2(2:N) = (Slhalf2(2:N) + Slhalf2(1:N-1))/2;

Slhalf2(isnan(Slhalf2) == 1) = 1;                                          % Setting NaN values to initial cross-sectional area
S_hat2(isnan(S_hat2) == 1) = 1;                                            % Setting NaN values to initial cross-sectional area

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


alpha0_p = (2*E*G0)/(1+(E*G0));                                            % Alpha0 loss coefficient for determining p
alphaq_p = (2*E*G1)/(1+(E*G0));                                            % Alpha q loss coefficient for determining p

xi_p = -(2*G1*k)/(2 + G0*k);                                               % vi loss coefficient for determining p0
eps_p = (2 - G0*k)/(2 + G0*k);                                             % epsilon loss coefficient for determining p0
nyu_p = (G0*k)/(2 + G0*k);                                                 % Nyu loss coefficient for determining p0

% ------------------------------------------------------------------------ %
%                              Lip Radiation
% ------------------------------------------------------------------------ %

lam2 = 2*lambda^2;                                                         % Lambda squared parameter

alpha0 = (L)/(0.8216*sqrt((S0*S1(N))/(pi)));                               % Loss parameter 2
alpha1 = 1/(2*0.8216^2*gamma);                                             % Loss parameter 1

% ------------------------------------------------------------------------- %
%                             Defining input
% ------------------------------------------------------------------------- %

if strcmp(input_type,'impulse')
    x = zeros(NF, 1);                                                      % Impulse input
    x(1) = 1;
    
elseif strcmp(input_type,'glottal')
    t = 0:1/Fs:Tf;                                                         % Time vector
    f_in = (depth/2)*cos(2*pi*fmod*t') + fmax-(depth/2);                    % frequency vector
    phase_in = cumsum(f_in/Fs);                                            % integral of frequency
    x = sin(2*pi*phase_in);                                                % generating pitch-varying sinusoid
    epsilon = 0.8;
    x = abs(x-epsilon)/1-epsilon;                                          % Applying abs function to get glottal waveform shape
    x(1:end-1) = x(1:end-1).*(1-breathiness) +...                          % Applying breathiness
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
D2(N+1,N) = 0;                                                             % Boundary condition
% ------------------------------------------------------------------------- %
%                             Initialising
% ------------------------------------------------------------------------- %

p = zeros(N+1, 1); p1 = zeros(N+1,1); p2 = zeros(N+1,1);                   % Initialising vectors for pressure

p0 = zeros(N+1,1); p01 = zeros(N+1,1); ptil = zeros(N+1,1); ...            % Initialising phat,p0
    ptil1 = zeros(N+1,1);

v = zeros(N,1); v1 = zeros(N,1); vpr = zeros(N,1); vpr1 = zeros(N,1);     % Initialising velocity vectors
v2 = zeros(N,1);

out = zeros(NF,1);                                                         % Initialing output vector

% ------------------------------------------------------------------------- %
%                             Update scheme
% ------------------------------------------------------------------------- %
tic;
for n= 1:NF
    
    S_hat = (1-interp(n))*S_hat1 + interp(n)*S_hat2;                       % Interpolate area functions for d_bar domain
    Slhalf = (1-interp(n))*Slhalf1 + interp(n)*Slhalf2;                    % Interpolate area functions for d domain 
    
    beta_p = (((k*rho*(gamma^2)))./((S_hat.*h)*((1+k*L))));                % Beta loss coefficient for pressure
    
    q1 = (alpha0*lambda^2*h*((3*S_hat(N) - S_hat(N-1))/(2*S_hat(N))));     % Boundary coefficient
    q2 = (alpha1*lambda^2*h*((3*S_hat(N) - S_hat(N-1))/(2*k*S_hat(N))));   % Boundary coefficient
    
    r1 = lam2 / (1+q1+q2);                                                 % defines radiation
    r2 = -(1+q1-q2)/(1+q1+q2);                                             % defines radiation
    
    v = alpha_v.*v - beta_v.*(D1*p)...
        + alpha_qv.*vpr;                                                   % Updating velocity
    
    vpr = tau_v*vpr + xi_v*(v + v1);                                       % Update v' 
    
    v(1) = x(n);                                                           % Excite system with input 
    
    p = alpha_p*p - ...
        beta_p.*(D2*(v.*Slhalf))...
        + alpha0_p*p0 + alphaq_p*ptil;                                     % Updating pressure
    
    p0 = eps_p*p01 + nyu_p*(p+p1)...                                       % Update p0
        + xi_p*ptil1;
    
    ptil = -ptil1 + ...
        (p + p1 - p0 - p01);                                               % Update pvhat
    
    p(N) = r1*p1(N-1) + r2*p2(N);
    
    if mod(n,anim) == 0                                                    % Drawing the pressure outputs
        if n == 1
            figure
            plt1 = plot(p(2:N));
            xlim([0 23])                                                   % Limiting axis
            ylim([-5*(10^(zoom_order)) 5*(10^(zoom_order))])
        else
            set(plt1, 'ydata', p)
            drawnow
        end
    end
    
    out(n) = p(N);                                                         % Reading output pressure difference
    v2 = v1; v1 = v; vpr1 = vpr; p01 = p0; ptil1 = ptil;p2=p1; p1 = p;
    
end
toc
% ------------------------------------------------------------------------- %
%                             Plotting outputs
% ------------------------------------------------------------------------- %


out = out/max(abs(out));                                                   % Normalise output



% Function to set the length of the vocal tract
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


% Function to give the area functions for two vowels
function [S1, S2] = getAreas(vowel)
if strcmp(vowel,'ao')
    
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('o');
elseif strcmp(vowel, 'ai')
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('i');
elseif strcmp(vowel,'ae')
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('e');
elseif strcmp(vowel,'au')
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('u');
elseif strcmp(vowel,'aa')
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('a');
elseif strcmp(vowel,'ea')
    S1 = getS_choudhury('e');
    S2 = getS_choudhury('a');
elseif strcmp(vowel, 'ee')
    S1 = getS_choudhury('e');
    S2 = getS_choudhury('e');
elseif strcmp(vowel,'ei')
    S1 = getS_choudhury('e');
    S2 = getS_choudhury('i');
elseif strcmp(vowel,'eo')
    S1 = getS_choudhury('e');
    S2 = getS_choudhury('o');
elseif strcmp(vowel,'eu')
    S1 = getS_choudhury('e');
    S2 = getS_choudhury('u');
elseif strcmp(vowel,'ia')
    S1 = getS_choudhury('i');
    S2 = getS_choudhury('a');
elseif strcmp(vowel, 'ie')
    S1 = getS_choudhury('i');
    S2 = getS_choudhury('e');
elseif strcmp(vowel,'ii')
    S1 = getS_choudhury('i');
    S2 = getS_choudhury('i');
elseif strcmp(vowel,'io')
    S1 = getS_choudhury('i');
    S2 = getS_choudhury('o');
elseif strcmp(vowel,'iu')
    S1 = getS_choudhury('i');
    S2 = getS_choudhury('u');
elseif strcmp(vowel,'oa')
    S1 = getS_choudhury('o');
    S2 = getS_choudhury('a');
elseif strcmp(vowel, 'oe')
    S1 = getS_choudhury('a');
    S2 = getS_choudhury('o');
elseif strcmp(vowel,'oi')
    S1 = getS_choudhury('o');
    S2 = getS_choudhury('i');
elseif strcmp(vowel,'oo')
    S1 = getS_choudhury('o');
    S2 = getS_choudhury('o');
elseif strcmp(vowel,'ou')
    S1 = getS_choudhury('o');
    S2 = getS_choudhury('u');
elseif strcmp(vowel,'ua')
    S1 = getS_choudhury('u');
    S2 = getS_choudhury('a');
elseif strcmp(vowel, 'ue')
    S1 = getS_choudhury('u');
    S2 = getS_choudhury('e');
elseif strcmp(vowel,'ui')
    S1 = getS_choudhury('u');
    S2 = getS_choudhury('i');
elseif strcmp(vowel,'uo')
    S1 = getS_choudhury('u');
    S2 = getS_choudhury('o');
elseif strcmp(vowel,'uu')
    S1 = getS_choudhury('u');
    S2 = getS_choudhury('u');
else
    error('input a combination of vowels please');
end
end

end
