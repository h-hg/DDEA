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
    %%%%    Authors:    Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty
    %%%%    University of Surrey, UK and Taiyuan University of Science and Technology, China.
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty, Offline data-driven evolutionary optimization using selective surrogate ensembles, IEEE Transactions on Evolutionary Computation, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
POP=lhsdesign(n,c).*(ones(n,1)*(bu-bd))+ones(n,1)*bd;

end