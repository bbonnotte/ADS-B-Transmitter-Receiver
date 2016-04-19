function [TEB ] = tauxDerreurBinaire(sb,b)
% Calcul du TEB
    % Input : sb, b
    % Outpout : TEB

N = length(sb);
    error = N - sum(b == sb);
    TEB = error / N;

end