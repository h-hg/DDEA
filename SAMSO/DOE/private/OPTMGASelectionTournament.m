function [parents, parentsScores] =OPTMGASelectionTournament(population,scores,NbTournaments)
%Function srgtsOPTMGASelectionTournament performs the GA tournament selection 
%operator. This operator performs a pair-wise ranking selection. Thus, for 
%example:
%
%     [PARENTS, PARENTSSCORES, NBPARENTS] = srgtsOPTMGASelectionTournament(OLDPOP,OLDSCORES,OPTIONS):
%     creates a new population with the selected individuals. OPTIONS is the
%     number of tournaments (which is equal the number of parents selected,
%     NBPARENTS).
%
%One can use this function out of the simpleOptimizer, even though it can not
%make any sense. Therefore, it is possible to manipulate the related option
%through the simpleOptimset, and then use this function inside "simpleOptimizer.m".
%
%Example: Create a SimpleToolBox option structure specifying the selection 
%function and also some related parameters.
%
%     options = simpleOptimset('GASelectionFcn', @srgtsOPTMGASelectionTournament, ...
%                              'GASelectionFcnOpts', 5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe Antonio Chegury Viana             Valder Steffen, Jr
% fchegury@yahoo.com                       vsteffen@mecanica.ufu.br
% http://fchegury.110mb.com
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PopSize   = length(scores);
tourns    = floor(rand(NbTournaments,PopSize)*PopSize) + 1; % Schedule of tournaments

% Determine the winner of the tournaments
[parentsScores idx] = min(reshape(scores(tourns),NbTournaments,PopSize));
parentsScores = parentsScores.';

% copy winners in to parents
parents = population(diag(tourns(idx,:)),:);

return
