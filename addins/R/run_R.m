function [slope1,slope2]=i_evqtlslopes(g012,expv)
    cd('C:\biodata\GEO\evQTL_Mapping\model3popgeno2');
    fid=fopen('r_fun/input.txt','w');
    for iix=1:length(g012)
        fprintf(fid,'%d\t%f\n',g012(iix),expv(iix));
    end
    fclose(fid);
    out=system('R CMD BATCH r_fun/runthisr.R');
    if out==0&&exist('r_fun/output.txt','file')
        z=importdata('r_fun/output.txt');
        slope1=z(1);
        slope2=z(2);
        delete('r_fun/output.txt');        
    else
        slope1=nan;
        slope2=nan;
    end    
    
    
    
