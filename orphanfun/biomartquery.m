function biomartquery()

x='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "HTML" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6"><Dataset name = "hsapiens_snp" interface = "default" ><Filter name = "refsnp" value = "rs12345"/><Attribute name = "refsnp_id" /><Attribute name = "chr_name" /><Attribute name = "chrom_start" /><Attribute name = "allele" /></Dataset></Query>';
urlFetch=sprintf('http://www.biomart.org/biomart/martservice?query=%s',...
    x);

    pagecontent=urlread(urlFetch);
