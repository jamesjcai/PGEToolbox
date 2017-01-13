function snp_vgrsidlist(rsidlist,popid)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

geno=[];
haplo=[];
for k=1:length(rsidlist)
    [geno1,mark1]=i_rsid2geno(rsidlist{k},popid);
    [geno2,mark2]=i_rsid2haplo(rsidlist{k},popid);
    geno=[geno,geno1];
    haplo=[haplo,geno2];
end
figure;
snp_vgview(geno);
figure;
snp_vhview(haplo);



function [geno1,mark1]=i_rsid2geno(marker,popid1)
       [s] = snp_downloadhapmap(marker,popid1);
       filename=tempname;
       [fid,Msg] = fopen(filename,'wt');
       if fid == -1, error(Msg); end
       fprintf(fid,'%s',s);
       fclose(fid);
       [geno1,mark1] = snp_readhapmap(filename);


function [geno1,mark1]=i_rsid2haplo(marker,popid1)
       [s] = snp_downloadhaplotype(marker,popid1);
       filename=tempname;
       [fid,Msg] = fopen(filename,'wt');
       if fid == -1, error(Msg); end
       fprintf(fid,'%s',s);
       fclose(fid);
       [geno1,mark1] = snp_readhaplotype(filename);

