%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ASSIGNMENT 6: PMMI
%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------- %
%                           Choir
% ------------------------------------------------------------------------- %

% This script consists of six calls to the time varying vocal tract that
% and sends it parameters depending upon the parameter choice in
% chords_Choudhury script. It is analogous to a choir of singers.



function tot = choir_Choudhury(opt, vib, chord, physical, sim)

% Checking length of inputs

if length(chord) == 1
    chord = repmat(chord,6);
elseif length(chord) ~= 6
    error('Please enter 6 values, one for each singer, or one value for all singers')
end

if length(physical.length) == 1
    physical.length = repmat(physical.length,6);
elseif length(physical.length) ~= 6
    error('Please enter 6 values, one for each singer, or one value for all singers')
end

if length(vib.depth) == 1
    vib.depth = repmat(vib.fdepth,6);
elseif length(vib.depth) ~= 6
    error('Please enter 6 values, one for each singer, or one value for all singers')
end

if length(vib.fmod) == 1
    vib.fmod = repmat(vib.fmod,6);
elseif length(vib.fmod) ~= 6
    error('Please enter 6 values, one for each singer, or one value for all singers')
end


%%%%% options
opts.anim = opt.animflag;                                                  % Animation flag
opts.animSpeed = opt.animspeed;                                            % Animation speed
opts.zoom = opt.animZoom;                                                  % Animation zoom
opts.type = sim.type;                                                      % Type of input waveform
opts.SR = sim.SR;                                                          % Sample rate
opts.Tf = sim.Tf;                                                          % Duration
opts.vowel = physical.vowel;                                               % Choice of vowel
opts.R0 = physical.R0;                                                     % Damping
opts.G0 = physical.G0;                                                     % Damping


%%%%% physical string parameters E
inp1.MIDI = chord(1);                                                      % Midi note
sim1.L = physical.length(1);                                               % length (m)
inp1.depth = vib.depth(1);                                                 % Depth of vibrato
inp1.fmod = vib.fmod(1);                                                   % Rate of vibrato
inp1.breath = physical.breathiness;                                        % Breathiness

%%%%% physical string parameters A 
inp2.MIDI = chord(2);
sim2.L = physical.length(2);                                               % length (m)
inp2.depth = vib.depth(2);
inp2.fmod = vib.fmod(2);
inp2.breath = physical.breathiness;

%%%%% physical string parameters D 
inp3.MIDI = chord(3);
sim3.L = physical.length(3);                                               % length (m)
inp3.depth = vib.depth(3);
inp3.fmod = vib.fmod(3);
inp3.breath = physical.breathiness;

%%%%% physical string parameters G 
inp4.MIDI = chord(4);
sim4.L = physical.length(4);                                               % length (m)
inp4.depth = vib.depth(4);
inp4.fmod = vib.fmod(4);
inp4.breath = physical.breathiness;

%%%%% physical string parameters B
inp5.MIDI = chord(5);
sim5.L = physical.length(5);                                               % length (m)
inp5.depth = vib.depth(5);
inp5.fmod = vib.fmod(5);
inp5.breath = physical.breathiness;

%%%%% physical string parameters E
inp6.MIDI = chord(6);
sim6.L = physical.length(6);                                               % length (m)
inp6.depth = vib.depth(6);
inp6.fmod = vib.fmod(6);
inp6.breath = physical.breathiness;

% ------------------------------------------------------------------------- %
%                             Outputting results
% ------------------------------------------------------------------------- %


sing1 = VT_TimeVarFunc_Choudhury(opts,sim1, inp1);
sing2 = VT_TimeVarFunc_Choudhury(opts,sim2, inp2);
sing3 = VT_TimeVarFunc_Choudhury(opts,sim3, inp3);
sing4 = VT_TimeVarFunc_Choudhury(opts,sim4, inp4);
sing5 = VT_TimeVarFunc_Choudhury(opts,sim5, inp5);
sing6 = VT_TimeVarFunc_Choudhury(opts,sim6, inp6);


A_pan = 1;                                                                 % Pan low E all the way to the left
B_pan = 4/5;                                                               % A slightly closer to right
C_pan = 3/5;                                                               % D fairly central panning
D_pan = 2/5;                                                               % Pan G fairly central, leaning right
E_pan = 1/5;                                                               % Pan B mainly on the left                                                              
F_pan = 0;                                                                 % Pan high E all the way to the right


tot(:,1) = A_pan*sing1 + B_pan*sing2 + C_pan*sing3 + ...
    D_pan*sing4 + C_pan*sing5 + F_pan*sing6;                               % Left channel

tot(:,2) = (1-A_pan)*sing1 + (1-B_pan)*sing2 + (1-C_pan)*sing3 + ...
    (1-D_pan)*sing4 + (1-E_pan)*sing5 + (1-F_pan)*sing6;                   % Right channel

end