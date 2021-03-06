% Author: Dr. Yuan SUN
% email address: yuan.sun@unimelb.edu.au OR suiyuanpku@gmail.com
% Modified by: Zhitao Huang
% email address: 2015160057@email.szu.edu.cn
%
% ------------
% Description:
% ------------
% This file is the entry point for the DECC-MDG algorithm on the CEC'2010 benchmark functions.
%


clear;
clc;

% set random seed
rand('state', sum(100*clock)); 
randn('state', sum(100*clock));
%warning('off' ,'Octave:divide-by-zero');

% problem dimension
D = 1000;

% population size
NP = 50;

% number of independent runs
runs = 25;

% number of fitness evaluations
Max_FEs = 3e6;

% for the benchmark functions initialization
global initial_flag;
problem=2010;
myfunc = 1:20;
addpath('benchmark2010');
addpath('benchmark2010/datafiles');
for func_num = myfunc 
    % load the FEs used by MDG in the decomposition process
    decResults = sprintf('./MergedDifferentialGrouping/results2010/F%02d',func_num);
    load (decResults);
    FEs = Max_FEs - FEs;
    
    % set the dimensionality and upper and lower bounds of the search space
    if(ismember(func_num, [1, 4, 7:9, 12:14, 17:20]))
        XRmin = -100*ones(NP,D); 
        XRmax = 100*ones(NP,D); 
        Lbound = XRmin;
        Ubound = XRmax;
    end
    if(ismember(func_num, [2, 5, 10, 15]))
        XRmin = -5*ones(NP,D); 
        XRmax = 5*ones(NP,D); 
        Lbound = XRmin;
        Ubound = XRmax;
    end
    if(ismember(func_num, [3, 6, 11, 16]))
        XRmin = -32*ones(NP,D); 
        XRmax = 32*ones(NP,D); 
        Lbound = XRmin;
        Ubound = XRmax;
    end
    
    Max_Gen = FEs/NP;

    VTRs = [];
    bestval = zeros(1,runs);
    for runindex = 1:runs
        % trace the fitness
        fprintf(1, 'Function %02d, Run %02d\n', func_num, runindex);
        filename = sprintf('trace2010/tracef%02d_%02d.txt',func_num, runindex);
        [fid, message] = fopen(filename, 'w');
        
        initial_flag = 0;
        % call the cmaescc algorithm
        [val]  = decc('benchmark_func', func_num, D, Lbound, Ubound, NP,Max_Gen, runindex,fid,problem);
        bestval(runindex) = val;
        fclose(fid);
    end
    
    meanval=mean(bestval);
    medianval=median(bestval);
    stdval=std(bestval);
    filename = sprintf('optimizationResults2010/f%02d.mat',func_num);
    save(filename,'bestval','meanval','medianval','stdval');

end

