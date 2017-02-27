function [pv,hetv2,fstv2,hetmd,fstmd]=fdist2test(genolist)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

geno=[];


nsubpop=length(genolist);
idvn=zeros(1,nsubpop);
for k=1:nsubpop
    g=genolist{k};
    idvn(k)=size(g,1);
    geno=[geno;g];
end

%[het,hetv]=snp_heterozygosity(geno);
%[fstv]=snp_fst3(genolist,'cross');
%fstv(fstv<0)=0;
%cd('C:\myprojects\pgetoolbox\addins\fdist2');

fprintf('Running FDIST2 DATACAL...\n');
[hetv2,fstv2]=i_datacalrun(geno,idvn);
% cmd=datacal.exe
% outfile=data_fst_outfile
fprintf('done.\n');

fprintf('Running FDIST2 simulations...');
%ndem=100;
ndem = nsubpop;
npop=nsubpop;
fstv1=fstv2(fstv2~=-100);
fstv1(fstv1<0)=0;

fst0=mean(fstv1);
nsmp=round(mean(idvn)*2);
muttype=0;
nrep=100000;
[hetmd,fstmd]=fdist2run(ndem,npop,fst0,nsmp,muttype,nrep);
% cmd=fdist2.exe
% outfile=fdist2out.dat
fprintf('done.\n');



fprintf('Running FDIST2 CPLOT...');
[x]=i_cplotrun;
% cmd=cplotcmd.exe
% outfile=cplotcmdout.dat
fprintf('done.\n');


fprintf('Running FDIST2 PV...');
[pv]=i_pvrun;
% cmd=pv.exe
% outfile=pvout.dat
% infile1=fdist2out.dat
% infile2=data_fst_outfile
fprintf('done.\n');


fprintf('Plotting FDIST2 Result...');

plot(hetmd,fstmd,'.','markersize',5,'color',[0.75294     0.75294     0.75294]);
hold on
plot(x(:,1),x(:,3),'-')
plot(x(:,1),x(:,[2 4]),'-b')
ylabel('F_{ST}')
xlabel('Het')

pv3=pv(fstv2~=-100);
hetv3=hetv2(fstv2~=-100);
fstv3=fstv2(fstv2~=-100);

plot(hetv3,fstv3,'.b')
plot(hetv3(pv3>0.95),fstv3(pv3>0.95),'.r')
plot(hetv3(pv3<0.05),fstv3(pv3<0.05),'.g')

xlim([0 0.52])
%for k=1:length(fstv2)
%text(hetv2(k),fstv2(k),sprintf('%d',k));
%end


hold off


fprintf('done.\n');





%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [x]=i_cplotrun
oldpath=pwd;
cdpge;
cd('addins/fdist2');

if ispc
    cplotcmd='cplotcmd.exe';
elseif ismac
    cplotcmd='./cplotcmd_mac';
else
    cplotcmd='./cplotcmd_unix';
end

try
[status] = system(cplotcmd);
catch
    warning('CPLOTCMD has not been compiled properly on this machine.');
    return;
end
if (status>0)
	error('FDIST2RUN: Runtime error.')
else
	%x=textread('cplotcmdout.dat');

        fid = fopen('cplotcmdout.dat');
        C = textscan(fid,'%f%f%f%f','Delimiter',' ','CollectOutput',1,...
            'TreatAsEmpty',{'1.#QNAN0'});
        fclose(fid);
        x=C{1};
        [a,b]=find(isnan(x));
        x(unique(a),:)=[];
end
cd(oldpath);




function [pv]=i_pvrun
oldpath=pwd;
cdpge;
cd('addins/fdist2');

if ispc
    pvcmd='pv.exe';
elseif ismac
    pvcmd='./pv_mac';
else
    pvcmd='./pv_unix';
end

[status] = system(pvcmd);
if (status>0)
	error('PV: Runtime error.')
else
	x=textread('pvout.dat');
    hetv=x(:,1);
    fstv=x(:,2);
    pv=x(:,4);
  %  hetv=hetv(fstv~=-100);
  %  fstv=fstv(fstv~=-100);  
end
cd(oldpath);




function [hetv,fstv]=i_datacalrun(geno,idvn)
oldpath=pwd;
cdpge;
cd('addins/fdist2');
subpopidx=[];
for k=1:length(idvn)
    subpopidx=[subpopidx;k*ones(idvn(k),1)];
end
fprintf('\tWriting to FDIST2 Input file...')
snp_writefdist2(geno,'infile',subpopidx);
fprintf('......done.')

[status] = system('datacal');
if (status>0)
	error('DATACAL: Runtime error.')
else
	x=textread('data_fst_outfile');
    hetv=x(:,1);
    fstv=x(:,2);
end
cd(oldpath);


