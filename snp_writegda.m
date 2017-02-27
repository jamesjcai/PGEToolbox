function [status] = snp_writegda(genodata,filename,subpopidx)
%snp_writegda - save genotype data for GDA
% [status] = snp_writegda(genodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

status=0;
[smpln,markn]=snp_samplen(genodata);

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
end
if ~(v==1), return; end


if nargin < 2
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), status=0; return; end
	filename=[pathname,filename];
end

fid1=fopen([filename,'.nex'],'w');

%pos=markinfo.pos;
%rsid=markinfo.rsid;

%First line: any character. Use this line to store information about your
%data.

fprintf(fid1, '#NEXUS\n\n');
fprintf(fid1, '[! \nSample data set on p. 338-339 of Weir, B. S. 1990. Genetic\n');
fprintf(fid1, 'Data Analysis. Sinauer, Sunderland, Mass.\n]\n\n');
fprintf(fid1, 'begin gdadata;\n');
fprintf(fid1, '\tdimensions nloci=%d npops=2;\n',markn);
fprintf(fid1, '\tformat tokens labels missing=? datapoint=standard;\n');

fprintf(fid1, '\tlocusallelelabels\n');
for k=1:markn
    fprintf(fid1, '\t\t%d ''locus-%d'',\n',k,k);
end
fprintf(fid1, '\t;\n\tMATRIX\n');

oldpop=0;
indc=1;
popnum=1;
for i=1:smpln
    if subpopidx(i)~=oldpop
        if popnum>1
            fprintf(fid1, '\t\t,\n');
        end
        fprintf(fid1, '\tPop_%d:\n',popnum);
        oldpop=subpopidx(i);
        indc=1;
        popnum=popnum+1;
    end
    fprintf(fid1, '\t\t_%d_ ',indc);
    indc=indc+1;
for j=1:2:markn*2
    if (genodata(i,j)==5 | genodata(i,j+1)==5)
        if j==markn*2-1
        fprintf(fid1, '?/?');            
            else
         fprintf(fid1, '?/? ');
        end
    else
        if j==markn*2-1
          fprintf(fid1, '%d/%d',genodata(i,j),genodata(i,j+1));
         else
        fprintf(fid1, '%d/%d ',genodata(i,j),genodata(i,j+1));
        end
    end    
end
    fprintf(fid1, '\n');
end
fprintf(fid1, '\t;\nEND;\n\n');


fprintf(fid1,'[\n');
fprintf(fid1,'Below is an example of a gda block. GDA blocks provide a means for\n');
fprintf(fid1,'running GDA in batch mode, specifying commands directly in the data\n');
fprintf(fid1,'file itself. This is a useful way to document your analyses. The\n');
fprintf(fid1,'GDA block below has been "turned off" by changing the name from "gda"\n');
fprintf(fid1,'to something the GDA program does not recognize (here, "_gda"). If\n');
fprintf(fid1,'GDA does not recognize a block name, it simply skips over it.\n');
fprintf(fid1,']\n');
fprintf(fid1,'begin _gda; [the underscore character prevents gda from recognizing this block]\n');
fprintf(fid1,'	log file=diploid.txt replace;\n');
fprintf(fid1,'	stats ppl al ap he ho samplesize f indiv ? est;\n');
fprintf(fid1,'	useloci 2-4;\n');
fprintf(fid1,'	fstats indivalleles ss ms vc coancestry noassumehw ? est;\n');
fprintf(fid1,'	bootloci rseed=1234567 nreps=9999 ci=95 ?;\n');
fprintf(fid1,'	jackpops;\n');
fprintf(fid1,'	jackloci;\n');
fprintf(fid1,'	gdist above=nei78iden below=nei78dist cluster=nj showmatrix drawphenogram nouselinedraw ? est;\n');
fprintf(fid1,'	usepop 2;\n');
fprintf(fid1,'	exact nruns=3200 subsets upto=3 noonlyhets measure=fisher missings=discard permute="ppppp" est;\n');
fprintf(fid1,'	log stop;\n');
fprintf(fid1,'end;\n');



fclose(fid1);
status=1;