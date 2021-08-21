% ------------
% Description:
% ------------
% This function is used to read the differential grouping datafiles and
% form the subcomponents which are used by the cooperative co-evolutionary
% framework. This function is called in decc.m.
%
%--------
% Inputs:
%--------
%    fun : the function id for which the corresponding grouping should
%          be loaded by reading the differential grouping datafiles.

function group = diff_grouping(fun,problem)   
    if(problem==2010)
        filename=sprintf('./MergedDifferentialGrouping/results2010/F%02d',fun);
    else
        filename=sprintf('./MergedDifferentialGrouping/results2013/F%02d',fun);
    end
    load(filename);
    group = nonseps;
    if(~isempty(seps))
        group = {group{1:end} seps};
    end
end
