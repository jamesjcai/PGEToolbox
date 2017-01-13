function [gwas]=gwascatalog(updatefile)

if nargin<1
	updatefile=false;
end

if ~updatefile
	load('gwascatalog','gwas');
else

	fullURL='http://www.genome.gov/admin/gwascatalog.txt';
	str=urlread(fullURL);
	temp=textscan(str,'%s','delimiter','\n');
	tags=textscan(temp{1}{1},'%s','delimiter','\t');
	content=textscan(str,repmat('%s',1,length(tags{1})),'delimiter','\t');

	disease={};
	gene={};
	allele={};

	for i=1:length(content{1})
	    tempgene=textscan(strrep(content{14}{i},' ',''),'%s','delimiter',',');
	    tempgene=tempgene{1};
	    disease=[disease; repmat(content{8}(i),length(tempgene),1)];
	    gene=[gene; tempgene];
	    allele=[allele; repmat(content{21}(i),length(tempgene),1)];
	end

	gwas.disease=disease;
	gwas.gene=gene;
	gwas.allele=allele;
	olddir=pwd;
        cdpge; 
	save('gwascatalog','gwas');
	cd(olddir);
end
