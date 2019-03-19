function [varargout] = size(X, DIM)
%size (overloaded)
%
% Author: Alexandre 
if nargin >= 2
	varargout{1} = size(X.data(1).value, DIM);
else
	%if nargout <= 1
		[varargout{1:nargout}] = size(X.data(1).value);
	%else
	%	[m, n] = size(X.data(1).value);
	%end
end
