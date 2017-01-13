% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


popcodes={'CEU','YRI','JPT+CHB'};

regionlist=[1   6095987	6096027
1	16509334	16509438
1	20792704	20792731
1	22273233	22273260];

genenamelist={'CHD5','FBXO42'};

k500=500000;

reg=[regionlist(:,1),...
regionlist(:,2)-500000,...
regionlist(:,3)+500000];

pk=1;
rk=1;

popcode=popcodes{pk};
query=sprintf('chr%d:%d-%d',reg(rk,:));
filename=sprintf('%s_%s_chr%d_%d_%d',popcode,genenamelist{rk},regionlist(rk,:));

%%
[s0] = snp_downloadhapmap(query,popcode);
fid=fopen(filename,'w');
fprintf(fid,s0);
fclose(fid)
