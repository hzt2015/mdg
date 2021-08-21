% Author: Zhitao Huang
% Email:  2015160057@email.szu.edu.cn 
%
% ------------
% Description:
% ------------
% This file is the entry point for running the merged differential 
% gropuing algorithm on the CEC'2013 benchmark functions.

clear all;

func = 1:15;
for func_num = func    
    t1 = [13 14];
    t2 = [1 4 7 8 11 12 15];
    t3 = [2 5 9];

    if (ismember(func_num, t1))
        D=905;
        lb = -100;
        ub = 100;
    elseif (ismember(func_num, t2))
        D=1000;
        lb = -100;
        ub = 100;
    elseif (ismember(func_num, t3))
        D=1000;
        lb = -5;
        ub = 5;
    else
        D=1000;
        lb = -32;
        ub = 32;
    end

    opts.lbound  = lb;
    opts.ubound  = ub;
    opts.dim     = D;
    %In order to prevent the special value of the upper and lower boundary 
    %from causing identification failure,the algorithm adjusts the varible value 
    %before perturbation to base, and the perturbation value is sigma
    opts.base    =(ub + lb)/2-(ub-lb)/8;
    opts.sigma=(ub-lb)/4;


    addpath('cec2013');
    addpath('cec2013/datafiles');
    global initial_flag;
    initial_flag = 0;

    [seps, nonseps, FEs] = MDG('benchmark_func', func_num, opts);
    filename = sprintf('./results2013/F%02d', func_num);
    save (filename, 'seps', 'nonseps', 'FEs', '-v7');
end

