function [status] = snp_writefdist2(genodata,filename,subpopidx)
%snp_writefdist2 - save genotype data for FDIST2
% [status] = snp_writefdist2(genodata,filename,subpopidx)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status=0;
[smpln,markn]=snp_samplen(genodata);
smpln=smpln/2;

if nargin<3
   subpopidx=ones(1,smpln);
   names={};
   for (k=1:smpln),
       names{k}=['Idv_',num2str(k)];
   end   
    [s,v] = choosebox('Name','Separate individuals','PromptString',...
        'Subpopulation 1:','SelectString','Subpopulation 2:',...
        'ListString',names); 
    subpopidx(s)=2;
if ~(v==1), return; end
end



if nargin < 2
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
end


%if strfind(filename,'.')>0
   fid1=fopen(filename,'w');
%else
%   fid1=fopen([filename,'.txt'],'w');
%end

%pos=markinfo.pos;
%rsid=markinfo.rsid;



%1/0 indicator. alleles by rows in data matrix (1), or pops by rows (0).
%No of pops
%no of loci

%no of alleles at locus 1
nsubpop=max(subpopidx(:));
fprintf(fid1, '0\n');
fprintf(fid1, '%d\n',nsubpop);
fprintf(fid1, '%d\n',markn);
fprintf(fid1, '\n');

for j=1:2:markn*2
    g=genodata(:,j:j+1);
    g2=g(:);     g2(g2==5)=[];
    ug2=unique(g2);
    [nalle]=length(ug2);   % no of alleles at locus J
    %fprintf(fid1, '%d\n',nalle);    
    %nalle=2;
    fprintf(fid1, '%d\n',nalle);
    
    for i=1:nsubpop
        gx=g(subpopidx==i,:);
        gx=gx(:); gx(gx==5)=[];
        for k=1:nalle
            fprintf(fid1, '%d\t', sum(gx==ug2(k)));
        end
        fprintf(fid1,'\n');
        %{        
        [minfreq,majfreq,minalle,majalle,mincount,majcount]=...
            snp_allefreq(g(subpopidx==i,:));
            if minalle<majalle
                   fprintf(fid1, '%d\t%d\n',mincount,majcount);
            else
                   fprintf(fid1, '%d\t%d\n',majcount,mincount);
            end
            %}
    end
    fprintf(fid1, '\n');    
end

fclose(fid1);
status=1;
