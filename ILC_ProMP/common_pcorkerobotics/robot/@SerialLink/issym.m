function res = issym(l)
%% ISSYM Checks wether Link or SerialLink object is a symbolic model.
%
%       res = issym(l)
%
%  Input::
%       l: Link or SerialLink class object (1x1)
%
%  Output::
%       res: true if l has symbolic variables as parameters.
%
%  Authors::
%        J�rn Malzahn   
%        2012 RST, Technische Universit�t Dortmund, Germany
%        http://www.rst.e-technik.tu-dortmund.de    
%

	res = issym(l.links(1));
	
end