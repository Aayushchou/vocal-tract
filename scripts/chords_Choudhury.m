clear all

% ------------------------------------------------------------------------- %
%                               Chord Sequencer 
% ------------------------------------------------------------------------- %

% This script calls the choir function with various parameters shaping the
% chords being played and how they are played. It calls the choir function
% 4 times to make a simple four chord bar of music. 

% ------------------------------------------------------------------------- %
%                          Global parameters
% ------------------------------------------------------------------------- %


opt.animspeed = 10;                                                        % Speed of animation 
opt.animZoom = 3;
opt.animflag = false;
sim.type = 'glottal';
sim.SR = 44100;
sim.Tf = 2;

% ------------------------------------------------------------------------- %
%                           First chord
% ------------------------------------------------------------------------- %

chord = [62 62 69 72 77 81];                                               % The MIDI notes being played
physical.breathiness = 0.1;                                                % Amount of breathiness
physical.length = [0.19 0.19 0.19 0.19 0.19 0.19];                         % Length of guitar: For actual size: use 0.648
physical.R0 = 1000;                                                        % Viscous losses
physical.G0 = 1000;                                                        % Thermal losses
physical.vowel = 'ai';                                                     % Choice of dipthong
vib.depth = [1, 1, 1, 1, 1, 1];                                            % Depth of vibrato 
vib.fmod  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];                                % Rate of vibrato
sim.SR = 44100;                                                            % Sample rate

if length(chord) ~= 6
    error('Enter 6 values for chord! one for each string')
end

p1 = choir_Choudhury(opt, vib, chord, physical, sim);                      % Calling choir
    

% ------------------------------------------------------------------------- %
%                          Second Chord 
% ------------------------------------------------------------------------- %

chord = [55 62 65 71 74 79];                                               % The notes being played
physical.breathiness = 0.1;
physical.length = [0.19 0.19 0.19 0.19 0.19 0.19];  
physical.R0 = 1000;
physical.G0 = 1000;
physical.vowel = 'oi';
vib.depth = [1, 1, 1, 1, 1, 1];
vib.fmod  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
sim.SR = 44100;

if length(chord) ~= 6
    error('Enter 6 values for chord! one for each string')
end

p2 = choir_Choudhury(opt, vib, chord, physical, sim);

% ------------------------------------------------------------------------- %
%                         Third chord
% ------------------------------------------------------------------------- %

chord = [60 64 67 71 76 79];                                               % The notes being played
physical.breathiness = 0.1;
physical.length = [0.19 0.19 0.19 0.19 0.19 0.19];
physical.R0 = 1000;
physical.G0 = 1000;
physical.vowel = 'ie';
vib.depth = [1, 1, 1, 1, 1, 1];
vib.fmod  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
sim.SR = 44100;

if length(chord) ~= 6
    error('Enter 6 values for chord! one for each string')
end

p3 = choir_Choudhury(opt, vib, chord, physical, sim);

% ------------------------------------------------------------------------- %
%                       Fourth Chord  
% ------------------------------------------------------------------------- %

chord = [61 64 67 71 76 79];                                               % The MIDI notes being played
physical.breathiness = 0.1;
physical.length = [0.19 0.19 0.19 0.19 0.19 0.19];
physical.R0 = 1000;                                                        % Damping 
physical.G0 = 1000;
physical.vowel = 'ia';
vib.depth = [1, 1, 1, 1, 1, 1];
vib.fmod  = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
sim.SR = 44100;

if length(chord) ~= 6
    error('Enter 6 values for chord! one for each string')
end

p4 = choir_Choudhury(opt, vib, chord, physical, sim);

out = [p1 ; p2 ; p3 ; p4];                                                 % Concatenate output 


out(:,1) = out(:,1)/max(abs(out(:,1)));                                    % Normalise left channel
out(:,2) = out(:,2)/max(abs(out(:,2)));                                    % Normalise right channel

soundsc(out, sim.SR)
audiowrite('chords.wav', out, sim.SR)
