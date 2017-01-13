function [idvlist,poplist,conlist]=smplpop_mapping_1kgenomes
%ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel
%http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/integrated_call_samples.20101123.ALL.panel
%ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel
olddir=pwd;
cdpge;
cd('addins')
cd('gzip')

fid=fopen('phase1_integrated_calls.20101123.ALL.panel','r');
[C]=textscan(fid,'%s%s%s%s%s%s');
fclose(fid);
    idvlist=C{1};
    poplist=C{2};
    conlist=C{3};
cd(olddir);


%{
fid=fopen('phase1_integrated_calls.20101123.ALL.panel','r');
[C]=textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines=C{1};

if length(lines)~=1092
    error('XXX');
end
idvlist=cell(1092,1);
poplist=cell(1092,1);
for k=1:length(lines)
    idvlist{k}=lines{k}(1:7);
    poplist{k}=lines{k}(9:11);
end
cd(olddir);
%}