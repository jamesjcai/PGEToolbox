function [fv,p1,p2,n1,n2]=snp_fstquery(query,popid1,popid2,ishapmap3)
% [fv,p1,p2,n1,n2]=snp_fstquery(query,popid1,popid2,ishapmap3)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2, popid1='CEU'; end
if nargin<3, popid2='YRI'; end
if nargin<4, ishapmap3=false; end

if (ishapmap3)
	downloadfun=@snp_downloadhapmap3;
else
	downloadfun=@snp_downloadhapmap;
end

p1=[]; p2=[]; n1=[]; n2=[];

fprintf('Downloading GENODATA from population 1 %s....\n',popid1);
       if (ishapmap3)
           [s,geno1,mark1]=downloadfun(query,popid1);
       else
           [geno1,mark1]=downloadfun(query,popid1);
       end
       if isempty(geno1)
           filename=tempname;
           [fid,Msg]=fopen(filename,'wt');
           if fid==-1, error(Msg); end
           fprintf(fid,'%s',s);
           fclose(fid);
           [geno1,mark1] = snp_readhapmap(filename);
       end
       if isempty(geno1), fv=[]; return; end
       if (~ishapmap3)
           [geno1]=snp_breaktrio(geno1,mark1);
       else
           %warning('Not yet SNP_BREAKTRIO')
       end

fprintf('Downloading GENODATA from population 2 %s....\n',popid2);
       if (ishapmap3)
           [s,geno2,mark2]=downloadfun(query,popid2);
       else
           [geno2,mark2]=downloadfun(query,popid2);
       end
       if isempty(geno2)   
           filename=tempname;
           [fid,Msg] = fopen(filename,'wt');
           if fid == -1, error(Msg); end
           fprintf(fid,'%s',s);
           fclose(fid);
           [geno2,mark2] = snp_readhapmap(filename);
       end
       if isempty(geno2), fv=[]; return; end
       if (~ishapmap3)
        [geno2]=snp_breaktrio(geno2,mark2);
       else
        %warning('Not yet SNP_BREAKTRIO')
       end


disp('Extracting common markers....')
x=intersect(mark1.pos,mark2.pos);
[a]=find(ismember(mark1.pos,x));
[b]=find(ismember(mark2.pos,x));
[geno1,mark1]=snp_pickmarker(geno1,mark1,a);
[geno2,mark2]=snp_pickmarker(geno2,mark2,b);


if (isempty(geno1) || isempty(geno2))
    fv=[];
    p1=[]; p2=[]; n1=[]; n2=[];
else
    [fv]=snp_fst(geno1,geno2);
    [p1,p2,n1,n2]=i_subpopallelefreq(geno1,geno2);
end



function [p1,p2,n1,n2]=i_subpopallelefreq(geno1,geno2)

s=2;                    % for SNP pairs, num of subpoulations, s=2
n1=sum(geno1(:)~=5);    % hiploid so n1 = n1*2;
n2=sum(geno2(:)~=5);
n=n1+n2;

%nc = (1/(s-1))*((n1+n2)-(n1^2+n2^2)/(n1+n2));
nc = (n-(n1^2+n2^2)/n);

[p1,alle1]=snp_maf(geno1);
[p2,alle2]=snp_maf(geno2);
if alle1~=alle2
    p2=1-p2;
end

