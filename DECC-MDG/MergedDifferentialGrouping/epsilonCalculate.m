function epsilon=epsilonCalculate(fp1,fp2,fp3,fp4,dim)
    FS = cat(1,abs(fp1),abs(fp2),abs(fp3),abs(fp4));
    Fmax=max(FS,[],1);     
    c=0.003;
    epsilon=c*realpow(2,-52)*Fmax*dim;
   
end