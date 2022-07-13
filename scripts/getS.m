% ------------------------------------------------------------------------ %
%                 Getter function: Cross-sectional area
% ------------------------------------------------------------------------ %

% This function inputs a character that indicates a vowel that the user
% wants to hear and outputs the corresponding cross-sectional area, by
% calling the data parser function, where the cross-sectional areas for
% different vowels are stored. The number argument corresponds to the
% column in which the cross-sectional area for the particular vowel is
% stored. 


function[S] = getS_choudhury(vowel)

if strcmp(vowel,'i')
    S = VowelDataParser_choudhury(2);                                      % store column 2 of S data
    S = S(1:end-2,:);                                                      % Reshape
    S(1:end,1) = S(1:end,1)/max(abs(S(1:end,1)));                          % Renormalise
    
elseif strcmp(vowel, 'e')
    S = VowelDataParser_choudhury(3);                                      % store column 3 of S data
    S = S(1:end-2,:);                                                      % Reshape
    S(1:end,1) = S(1:end,1)/max(abs(S(1:end,1)));                          % Renormalise
    
elseif strcmp(vowel, 'a')
    S = VowelDataParser_choudhury(7);                                      % store column 7 of S data
    
elseif strcmp(vowel, 'o')
    S = VowelDataParser_choudhury(9);                                      % store column 9 of S data
        
elseif strcmp(vowel, 'u')
    S = VowelDataParser_choudhury(12);                                     % store column 11 of S data

else
    error('Please pick a one of the available vowels (a,e)')
end

end
