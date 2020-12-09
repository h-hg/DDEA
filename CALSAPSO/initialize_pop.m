function [ POP ] = initialize_pop(n,c,bu,bd)
% Usage: [ POP ] = initialize_pop(n,c,bu,bd)
%
% Input:
% bu            -Upper Bound
% bd            -Lower Bound
% c             -No. of Decision Variables
% n             -Population Scale
%
% Output: 
% POP           -Initial Population
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, John Doherty
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:
%Handing Wang, Yaochu Jin, John Doherty, Committee-based Active Learning for Surrogate-Assisted Particle Swarm Optimization of Expensive Problems, IEEE Transactions on Cybernetics, vol.47, no.9, pp.2664-2677, 2017.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
POP=lhsdesign(n,c).*(ones(n,1)*(bu-bd))+ones(n,1)*bd;

end