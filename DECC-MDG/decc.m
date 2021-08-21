% Author: Dr. Zhyenu Yang
% Modified by: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
%
% ------------
% Description:
% ------------
% This file is an implementation of cooperative co-evolution which
% uses SaNSDE algorithm as subcomponent optimizer.
%
% -----------
% References:
% -----------
% Omidvar, M.N.; Li, X.; Mei, Y.; Yao, X., "Cooperative Co-evolution with
% Differential Grouping for Large Scale Optimization," Evolutionary Computation,
% IEEE Transactions on, vol.PP, no.99, pp.1,1, 0
% http://dx.doi.org/10.1109/TEVC.2013.2281543
%
% --------
% License:
% --------
% This program is to be used under the terms of the GNU General Public License 
% (http://www.gnu.org/copyleft/gpl.html).
% Author: Mohammad Nabi Omidvar
% e-mail: mn.omidvar AT gmail.com
% Copyright notice: (c) 2013 Mohammad Nabi Omidvar


function [bestval] = decc(fname, func_num, dim, Lbound, Ubound, popsize, itermax, runindex, fid,problem)

    % for fitness trace
    tracerst = [];

    % the initial population
    pop = Lbound + rand(popsize, dim) .* (Ubound-Lbound);

    val = feval(fname, pop, func_num);
    [bestval, ibest] = min(val);%Return the minimum value of each column and the row number corresponding to the value (starting from 1)
    bestmem = pop(ibest, :);%Get the best fit individual


    % the initial crossover rate for SaNSDE
    group = {};
    ccm = 0.5;%CRm交叉率的高斯分布随机函数的平均数，每25代更新一次
    sansde_iter = 1;
    Cycle = 0;
    iter = 0;
    delta = 0;

    FE = 0;
    
    group = diff_grouping(func_num,problem);%返回差分分组的分组情况
    group_num = size(group, 2);%返回分组数
        
    
    display = 0;
    frequency = 100;

    while (iter < itermax)

        for i = 1:group_num
            oneitermax = sansde_iter;
            if (iter + oneitermax >= itermax)
                oneitermax = itermax - iter;
            end
            if (oneitermax == 0)
                break;
            end

            dim_index = group{i};%返回第i个子组件的维度（决策变量）索引
            subpop = pop(:, dim_index); 
            subLbound = Lbound(:, dim_index);        
            subUbound = Ubound(:, dim_index);

                [subpopnew, bestmemnew, bestvalnew, tracerst, ccm] = sansde(fname, func_num, dim_index, subpop, bestmem, bestval, subLbound, subUbound, oneitermax, ccm,Cycle);

                FE = FE + popsize;

                iter = iter + oneitermax;
                
                %输出格式为实数，科学计算法形式
                fprintf(fid, '%e\n', tracerst); 

                
                pop(:, dim_index) = subpopnew;
                bestmem = bestmemnew;
                bestval = bestvalnew;

                nonsep = [];
                sep = [];

                if(display == 1)
                    fprintf(1, 'Cycle = %d, bestval = %e, Group = %d *\n',  Cycle, bestval, i);
                end

                if(iter > itermax)
                    break;
                end


        end

     
        val = feval(fname, pop, func_num);
        [best, ibest] = min(val);
        if (best < bestval)
            bestval = best;
            bestmem = pop(ibest, :);
        end

        Cycle = Cycle + 1;
        
        
    end
 
end
