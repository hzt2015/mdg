% Author: Zhitao Huang
% Email:  2015160057@email.szu.edu.cn 
%
% ------------
% Description:
% ------------
% This file is the entry point for running the merged differential 
% gropuing algorithm on the CEC'2010 benchmark functions.


clear all;
func=1:20;

for func_num=func
    t1 = [1 4 7 8 9 12 13 14 17 18 19 20];
    t2 = [2 5 10 15];
    t3 = [3 6 11 16];
    
    if (ismember(func_num, t1))
        lb = -100;
        ub = 100;
    elseif (ismember(func_num, t2))
        lb = -5;
        ub = 5;
    elseif (ismember(func_num, t3))
        lb = -32;
        ub = 32;
    end

    opts.lbound  = lb;
    opts.ubound  = ub;
    opts.dim     = 1000;
    %In order to prevent the special value of the upper and lower boundary 
    %from causing identification failure,the algorithm adjusts the varible value 
    %before perturbation to base, and the perturbation value is sigma
    opts.base    = (ub + lb)/2-(ub-lb)/8;
    opts.sigma   = (ub-lb)/4;
    
    addpath('cec2010');
    addpath('cec2010/datafiles');
    global initial_flag;
    initial_flag = 0;
    
    [seps, nonseps, FEs] = MDG('benchmark_func', func_num, opts);
    filename = sprintf('./results2010/F%02d', func_num);
    save (filename, 'seps', 'nonseps', 'FEs', '-v7'); 
       
end    