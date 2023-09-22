function [hap, mrk] = i_hap_pickmarker(hap, mrk, idx)
if isempty(idx), return; end
hap = hap(:, idx);
mrk.rsid = mrk.rsid(idx);
mrk.pos = mrk.pos(idx);
mrk.maf = mrk.maf(idx);