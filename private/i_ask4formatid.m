function [id] = i_ask4formatid()

	ButtonName=questdlg('What kind of file format?', ...
			    'Select sequence format', ...
			    'FASTAS','PHYLIP','MAT','FASTAS');
	switch ButtonName,
	    case 'FASTAS', 
		id=1;          
	    case 'PHYLIP'
		id=2;    
	    case 'MAT'
		id=3;        
	    otherwise
		id=4;
	end