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
%    dims       : Variable vector currently to be grouped
%
%    fp1        : Function value without perturbation of all variables
%
%    perturbed_values   : Function value after perturbing each variable
%
%
% --------
% Outputs:
% --------
%    groups   : The grouping result of the current vector dims.
%    groups_perturb: Function value after perturbing dims.
%    FEs      : the number of fitness evaluations used.
%
% ------
function [groups,groups_perturb,FEs]=mergeGroup(fun,fun_number,options,dims,fp1,perturbed_values)
    dim=options.dim;
    base=options.base;
    sigma=options.sigma;
    FEs=0;
    dim_len=length(dims);
    groups={};
    
    %% Only one variable, or two variables
    if(dim_len==1 || dim_len==2) 
        %A variable, return directly
        if(dim_len==1)
            groups_perturb=perturbed_values(dims);
            groups{1}=dims;
        else    %Two variables, determine whether to interact
            p=base * ones(1,dim);
            p(dims)=base+sigma;
            fp=feval(fun,p,fun_number);
            FEs=FEs+1;
            
            delta1=perturbed_values(dims(1))-fp1;
            delta2=fp-perturbed_values(dims(2));
            epsilon=epsilonCalculate(fp1,perturbed_values(dims(1)),perturbed_values(dims(2)),fp,dim);
            
            groups_perturb=fp;
            if (abs(delta1-delta2)<epsilon)
                groups{1}=dims(1);
                groups{2}=dims(2);
            else
                groups{1}=dims';
            end
        end
        return;
    end
    
    %% The variable set dims is divided into two subsets
    median=floor(dim_len/2);
    Ldims=dims(1:median);
    Rdims=dims(median+1:dim_len);
    
    %Two subsets are grouped separately
    [Lgroups,Lgroups_perturb,LFEs]=mergeGroup(fun,fun_number,options,Ldims,fp1,perturbed_values);
    [Rgroups,Rgroups_perturb,RFEs]=mergeGroup(fun,fun_number,options,Rdims,fp1,perturbed_values);
    
    FEs=LFEs+RFEs;
    L_gnum=size(Lgroups,2);
    R_gnum=size(Rgroups,2);
    
    %% Merge non-separable subsets between two subset groups
    %Determine whether there is an interaction between two subset groups
    Lfp=Lgroups_perturb;
    Rfp=Rgroups_perturb;
    
    p=base * ones(1,dim);
    p(dims)=base+sigma;
    fp=feval(fun,p,fun_number);
    FEs=FEs+1;
    
    delta1=Lfp-fp1;
    delta2=fp-Rfp;
    epsilon=epsilonCalculate(fp1,Lfp,Rfp,fp,dim);
    groups_perturb=fp;
    
    
    %There are interactions between subset groups  
    if(abs(delta1-delta2)>epsilon) 
        if(L_gnum==1&&R_gnum==1)
            Rgroups{1}=[Lgroups{1} Rgroups{1}];
            groups=Rgroups;
            return ; 
        end
        
        %First find out the subset of the left subset group that interacts with the right subset group
        [Lgroup,LgroupIndexs,~,FE]=biSearch(fun,fun_number,options,Rdims',fp1,Rfp,fp,Lgroups,Lfp);
        FEs=FEs+FE;
        if(~isempty(LgroupIndexs))
            lnum=length(LgroupIndexs);
            %There is only a subset on the right
            if(R_gnum==1)
                Rgroups{1}=[cell2mat(Lgroup) Rgroups{1}];
                LgroupIndexs=sort(LgroupIndexs,'descend');
                for i=LgroupIndexs
                    Lgroups(i)=[];
                end
                groups=[Lgroups Rgroups];
                return ; 
            end
            Rinteract_Indexs=cell(1,lnum);
            Rperts=Rfp;%Perturbation function value of the variable set corresponding to the binary search tree on the right
            for i=1:lnum
                if(L_gnum==1)
                    fp_i=Lfp;
                    fp_iR=fp;
                else
                    p_i=base * ones(1,dim);
                    p_i(Lgroup{i})=base+sigma;
                    fp_i=feval(fun,p_i,fun_number);
                    p_iR=p_i; 
                    iR=[Lgroup{i}, Rdims'];
                    p_iR(iR)=base+sigma;
                    fp_iR=feval(fun,p_iR,fun_number);
                    FEs=FEs+2;
                end            
                [~,RgroupIndexs,Rperts,FE]=biSearch(fun,fun_number,options,Lgroup{i},fp1,fp_i,fp_iR,Rgroups,Rperts);
                FEs=FEs+FE;
                Rinteract_Indexs{i}=RgroupIndexs;
            end
            
            %Merge interaction group
            [merged_groups]=mergeInteractionGroup(LgroupIndexs,Rinteract_Indexs,Lgroups,Rgroups);
            LgroupIndexs=sort(LgroupIndexs,'descend');
            for j=LgroupIndexs
                Lgroups(j)=[];
            end
            list_r=[];
            for k=1:lnum
                list_r=union(list_r,Rinteract_Indexs{k});
            end
            list_r=sort(list_r,'descend');
            for j=list_r
                Rgroups(j)=[];
            end
            Rgroups=[Rgroups merged_groups];
        end
    end
        
    groups=[Lgroups Rgroups];    
           
end