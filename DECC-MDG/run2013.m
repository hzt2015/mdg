% Author: Dr. Yuan SUN
% email address: yuan.sun@unimelb.edu.au OR suiyuanpku@gmail.com
% Modified by: Zhitao Huang
% email address: 2015160057@email.szu.edu.cn
%
% ------------
% Description:
% ------------
% This file is the entry point for the DECC-MDG algorithm on the CEC'2013 benchmark functions.
%


clear;
% set random seed
rand('state', sum(100*clock)); 
randn('state', sum(100*clock));
%warning('off' ,'Octave:divide-by-zero');

% number of independent runs
runs = 25;

% population size
NP = 50;

% number of fitness evaluations
Max_FEs = 3e6;

% for the benchmark functions initialization
global initial_flag;
myfunc = 1:15;
addpath('benchmark2013');
addpath('benchmark2013/datafiles');
problem=2013;
for func_num = myfunc 
    % load the FEs used by MDG in the decomposition process
    decResults = sprintf('./MergedDifferentialGrouping/results2013/F%02d',func_num);
    load (decResults);
    FEs = Max_FEs - FEs;
    
    % set the dimensionality and upper and lower bounds of the search space
    if (ismember(func_num, [13,14]))
        D = 905;
        Lbound = -100.*ones(NP,D);
        Ubound = 100.*ones(NP,D);
    elseif (ismember(func_num, [1,4,7,8,11,12,15]))
        D = 1000;
        Lbound = -100.*ones(NP,D);
        Ubound = 100.*ones(NP,D);
    elseif (ismember(func_num, [2,5,9]))
        D=1000;
        Lbound = -5.*ones(NP,D);
        Ubound = 5.*ones(NP,D);
    else 
        D=1000;
        Lbound = -32.*ones(NP,D);
        Ubound = 32.*ones(NP,D);
    end
    
    Max_Gen = FEs/NP;

    VTRs = [];
    bestval = zeros(1,runs);
    for runindex = 1:runs
        % trace the fitness
        fprintf(1, 'Function %02d, Run %02d\n', func_num, runindex);
        filename = sprintf('trace2013/tracef%02d_%02d.txt',func_num, runindex);
        [fid, message] = fopen(filename, 'w');
        
        initial_flag = 0;
        % call the decc algorithm
        [val]  = decc('benchmark_func', func_num, D, Lbound, Ubound, NP,Max_Gen, runindex,fid,problem);
        bestval(runindex) = val;
        fclose(fid);
    end
    
    meanval=mean(bestval);
    medianval=median(bestval);
    stdval=std(bestval);
    filename = sprintf('optimizationResults2013/f%02d.mat',func_num);
    save(filename,'bestval','meanval','medianval','stdval');

end

