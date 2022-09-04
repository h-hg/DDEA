classdef INDIVIDUAL < handle
%INDIVIDUAL - The class of an individual.%%%%11111
%
%   This is the class of an individual. An object of INDIVIDUAL stores all
%   the properties including decision variables, objective values,
%   constraint violations, and additional properties of an individual.
%
% INDIVIDUAL properties:
%   dec         <read-only>     decision variables of the individual
%   obj         <read-only>     objective values of the individual
%   con         <read-only>     constraint violations of the individual
%   add         <read-only>     additional properties of the individual
%
% INDIVIDUAL methods:
%   INDIVIDUAL	<public>        the constructor, all the properties will be
%                               set when the object is creating
%   decs        <public>      	get the matrix of decision variables of the
%                               population
%   objs        <public>        get the matrix of objective values of the
%                               population
%   cons        <public>        get the matrix of constraint violations of
%                               the population
%   adds        <public>        get the matrix of additional properties of
%                               the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        dec;        % Decision variables of the individual
        obj;        % Objective values of the individual
        con;        % Constraint violations of the individual
        add;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = INDIVIDUAL(DEC,OBJ)
                 obj.dec=DEC;
                 obj.obj=OBJ;

        end
               %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.dec);
        end
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.obj);
        end
    end
end