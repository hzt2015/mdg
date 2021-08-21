% -------
% Inputs:
% -------
%    fun        : the function suite for which the interaction structure
%                 is going to be identified in this case benchmark_func
%                 of cec'2010 or cec'2013.
%    fun_number : the function number.
%    options    : this variable contains the options such as problem
%                 dimensionality, upper and lower bounds, and 
%                 parameters base and sigma used by merged differential grouping.   
%    left_group : The current subset of the interaction to be identified in Lgroups
%    Rgroups    : The right subset group
%    fp1        : Function value without perturbation of all variables
%    fp2        : Function value perturbing left_group.
%    fp_iR      : Function value perturbing left_group and Rgroups.
%    Rperts   : Perturbation function value of the variable set corresponding to the binary search tree on the right
%
% --------
% Outputs:
% --------
%    group   : The subset in Rgroups interacting with the current left subset left_group.
%    groupIndexs: Corresponding Index of Subset of Interaction.
%    Rperts   : Perturbation function value of the variable set corresponding to the binary search tree on the right
%    FEs      : the number of fitness evaluations used.
%
% ------

function [group,groupIndexs,Rperts,FEs]=biSearch(fun,fun_number,options,left_group,fp1,fp2,fp_iR,Rgroups,Rperts)
    base=options.base;
    sigma=options.sigma;
    dim=options.dim;
    FEs=0;
    R_gnum=size(Rgroups,2); % The number of subsets in Rgroups
    Rfp=Rperts(1);          % Perturbation function value of all variables in Rgroups
    groupsQueue={};         % Save subset groups of the breadth-first traversal strategy
    dataQueue={};           % Save the function value of the subset groups of the breadth-first traversal strategy
    groupIndexsQueue={};    % Save the index position of the subset in each node in Rgroups
    pQueue={};              % Save the basic decision vector value of the current node
    
    data=[fp1,fp2,Rfp,fp_iR];
    groupsQueue={Rgroups};
    dataQueue{1}=data;
    groupIndexsQueue{1}=1:R_gnum;
    pQueue{1}=base * ones(1,dim);
    node_orders=1;
    group={};
    groupIndexs=[];
    
    while(~isempty(groupsQueue))
        %Take out the current team leader
        cur_groups=groupsQueue{1};
        cur_data=dataQueue{1};
        cur_groupIndexs=groupIndexsQueue{1};
        cur_p=pQueue{1};
        cur_order=node_orders(1);
        groupsQueue(1)=[];
        dataQueue(1)=[];
        pQueue(1)=[];
        groupIndexsQueue(1)=[];
        node_orders(1)=[];
        delta1=cur_data(2)-cur_data(1);
        delta2=cur_data(4)-cur_data(3);
        epsilon=epsilonCalculate(cur_data(1),cur_data(2),cur_data(3),cur_data(4),dim);
%         There is an interaction between the current subset group and left_group, 
%         and the current subset group is divided into two small subset groups
        if(abs(delta1-delta2)>epsilon)
            cur_gnum=size(cur_groups,2);
            if(cur_gnum==1) 
                group={group{1:end} cur_groups{1}};
                groupIndexs(end+1)=cur_groupIndexs(1);
            else
                median=floor(cur_gnum/2);
                groups1=cur_groups(1:median);
                groupIndexs1=cur_groupIndexs(1:median);
                p_1=cur_p;
                p_1(cell2mat(groups1))=base+sigma;
                if(cur_order*2<=length(Rperts)&&Rperts(cur_order*2)~=0)
                    fp_1=Rperts(cur_order*2);
                else
                    fp_1=feval(fun,p_1,fun_number);
                    Rperts(cur_order*2)=fp_1;
                    FEs=FEs+1;
                end
                p_i1=p_1;
                p_i1(left_group)=base+sigma;
                fp_i1=feval(fun,p_i1,fun_number);
                FEs=FEs+1;
                
                data1=[cur_data(1),cur_data(2),fp_1,fp_i1];
                groupsQueue{end+1}=groups1; 
                dataQueue{end+1}=data1;
                groupIndexsQueue{end+1}=groupIndexs1;
                pQueue{end+1}=cur_p;
                node_orders(end+1)=cur_order*2;

                groups2=cur_groups(median+1:cur_gnum);
                groupIndexs2=cur_groupIndexs(median+1:cur_gnum);
                data2=[fp_1,fp_i1,cur_data(3),cur_data(4)];
                groupsQueue{end+1}=groups2;
                dataQueue{end+1}=data2;
                groupIndexsQueue{end+1}=groupIndexs2;
                pQueue{end+1}=p_1; 
                node_orders(end+1)=cur_order*2+1;
            end           
        end
    end
     
end