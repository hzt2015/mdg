function [merged_groups]=mergeInteractionGroup(LgroupIndexs,Rinteract_Indexs,Lgroups,Rgroups)
    merged_groups={};
    while(~isempty(LgroupIndexs))
        cur_index=LgroupIndexs(1);
        LgroupIndexs(1)=[];
        cur_interact_indexs=Rinteract_Indexs{1};
        Rinteract_Indexs(1)=[];
        
        groups=Lgroups{cur_index};
        len_l=length(LgroupIndexs);
        i=len_l;
        while i>0
            if(~isempty(intersect(cur_interact_indexs,Rinteract_Indexs{i})))
                cur_interact_indexs=union(cur_interact_indexs,Rinteract_Indexs{i});
                groups=[groups Lgroups{LgroupIndexs(i)}];
                Rinteract_Indexs(i)=[];
                LgroupIndexs(i)=[];
                len_l=len_l-1;
                i=len_l;
            else
                i=i-1; 
            end
        end
        for i=cur_interact_indexs
            groups=[groups Rgroups{i}];
        end
        merged_groups=[merged_groups groups];
    end
end