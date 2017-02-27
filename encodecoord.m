function [P,B] = encodecoord(selection)
%ENCODECOORD - genomic coordinates of ENCODE regions

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<1), selection=44; end

switch (selection)
    case (44)
P=[23 122609996 123109995; 23 152767492 154063081; 1 149424685 149924684; 2 234156564 234656627; 2 219985590 220485589; 2 51512209 52012208; 2 118011044 118511043; 4 118466104 118966103; 5 141880151 142380150; 5 131284314 132284313; 5 55871007 56371006; 6 132218540 132718539; 6 73789953 74289952; 6 108371397 108871396; 6 41405895 41905894; 7 89621625 90736048; 7 115597757 117475182; 7 26924046 27424045; 7 113720369 114720368; 7 125865892 127029088; 8 118882221 119382220; 9 130725123 131225122; 10 55153819 55653818; 11 130604798 131104797; 11 63940889 64440888; 11 4730996 5732587; 11 1699992 2306039; 11 115962316 116462315; 12 38626477 39126476; 13 29418016 29918015; 13 112338065 112838064; 14 52947076 53447075; 14 98458224 98958223; 15 41520089 42020088; 16 1 500000; 16 60833950 61333949; 16 25780428 26280428; 18 23719232 24219231; 18 59412301 59912300; 19 59023585 60024460; 20 33304929 33804928; 21 39244467 39744466; 21 32668237 34364221; 22 30133954 31833953];
B={
'ENr324 ','X','Random Picks';
'ENm006 ','X','Manual Picks:ChrX';
'ENr231 ','1','Random Picks';
'ENr131 ','2','Random Picks';
'ENr331 ','2','Random Picks';
'ENr112 ','2','Random Picks';
'ENr121 ','2','Random Picks';
'ENr113 ','4','Random Picks';
'ENr212 ','5','Random Picks';
'ENm002 ','5','Manual Picks:Interleukin';
'ENr221 ','5','Random Picks';
'ENr222 ','6','Random Picks';
'ENr223 ','6','Random Picks';
'ENr323 ','6','Random Picks';
'ENr334 ','6','Random Picks';
'ENm013 ','7','Manual Picks';
'ENm001 ','7','Manual Picks:CFTR';
'ENm010 ','7','Manual Picks:HOXA';
'ENm012 ','7','Manual Picks:FOXP2';
'ENm014 ','7','Manual Picks';
'ENr321 ','8','Random Picks';
'ENr232 ','9','Random Picks';
'ENr114 ','10','Random Picks';
'ENr312 ','11','Random Picks';
'ENr332 ','11','Random Picks';
'ENm009 ','11','Manual Picks:Beta';
'ENm011 ','11','Manual Picks:1GF2/H19';
'ENm003 ','11','Manual Picks:Apo';
'ENr123 ','12','Random Picks';
'ENr111 ','13','Random Picks';
'ENr132 ','13','Random Picks';
'ENr311 ','14','Random Picks';
'ENr322 ','14','Random Picks';
'ENr233 ','15','Random Picks';
'ENm008 ','16','Manual Picks:Alpha';
'ENr313 ','16','Random Picks';
'ENr211 ','16','Random Picks';
'ENr213 ','18','Random Picks';
'ENr122 ','18','Random Picks';
'ENm007 ','19','Manual Picks:Chr19';
'ENr333 ','20','Random Picks';
'ENr133 ','21','Random Picks';
'ENm005 ','21','Manual Picks:Chr21';
'ENm004 ','22','Manual Picks:Chr22'};

    case (10)

    % http://www.hapmap.org/downloads/encode1.html.en

           P = [2 51633239 52133238;
		2 234778639 235278638;
		4 118705475 119205474;
		7 26699793 27199792;
		7 89395718 89895717;
		7 126135436 126632577;
		8 118769628 119269627;
		9 127061347 127561346;
		12 38626477 39126476;
		18 23717221 24217220];
B = {'2p16.3','ENr112';'2q37.1','ENr131';'4q26','ENr113';'7p15.2','ENm010';'7q21.13','ENm013';'7q31.33','ENm014';'8q24.11','ENr321';'9q34.11','ENr232';'12q12','ENr123';'18q12.1','ENr213'};

    otherwise
	 error('Wrong Selection.')
end




