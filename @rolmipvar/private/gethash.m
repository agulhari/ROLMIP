function [index] = gethash(exponent, exptable, jump)
	% defined for convenience if gethash.c was not compiled into a mex file
	index = 0;
	for contsimplex = 1:(length(exponent)) 
		if (length(exptable{contsimplex}) > 0)
			[i, j] = find(exptable{contsimplex} - repmat(exponent{contsimplex},size(exptable{contsimplex},1), 1) == 0);
			index = index + (find(histc(i, 1:size(exptable{contsimplex}, 1)) == size(exptable{contsimplex},2)) - 1)*jump(contsimplex);
		end
	end
	index = index + 1;
end