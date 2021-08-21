% -------
% Inputs:
% -------
%    fun        : the function suite for which the interaction structure
%                 is going to be identified in this case benchmark_func
%                 of cec'2010 or cec'2013.
%
%    fun_number : the function number.
%
%    options    : this variable contains the options such as problem
%                 dimensionality, upper and lower bounds, and 
%                 parameters base and sigma used by merged differential grouping.
%
% --------
% Outputs:
% --------
%    sep      : a vector of all separable variables.
%    allgroups: a cell array containing all non-separable groups.
%    FEs      : the total number of fitness evaluations used.
%
% ------


function [seps, allgroups, FEs] = MDG(fun, fun_number, options)
    dim       = options.dim;
    base      = options.base;
    sigma     = options.sigma;
    dims      = [];%Save non-separable variables
    seps      = [];%Save separable variables
    allgroups = {};%save non-separable subcomponents
    FEs       = 0;
    perturbed_values=[];%Function value after perturbing each variable

    %% Find all the separable variables
    p1 = base * ones(1,dim);
    fp1 = feval(fun, p1, fun_number);
    p4 = (base+sigma) * ones(1,dim);
    fp4 = feval(fun, p4,fun_number);

    FEs = FEs + 2;

    for i =1:dim
        p2=p1;
        p2(i)=(base+sigma);
        fp2=feval(fun, p2, fun_number);
        
        p3=p4;
        p3(i)=base;
        fp3=feval(fun, p3, fun_number);
        
        FEs = FEs +2;

        delta1=fp2-fp1;
        delta2=fp4-fp3;
        perturbed_values=[perturbed_values;fp2];

        
%         adaptive epsilon
        epsilon=epsilonCalculate(fp1,fp2,fp3,fp4,dim); 
        if(abs(delta1-delta2) < epsilon)            
            seps = [seps ;i];
        else
            dims=[dims;i];
        end
    end
    
   %% grouping nonseparable variables
     if(length(dims)>1)
        [groups,~,FE]=mergeGroup(fun,fun_number,options,dims,fp1,perturbed_values);
        FEs = FEs + FE;
        %Find groups with only one variable, and merge the variables into seps vector
        for i=size(groups,2):-1:1
            if(length(groups{i})==1)
                seps = [seps ;groups{i}];
                groups(i)=[];
            end
        end
        allgroups=groups;
     end
end
