function [sys] = psys(this, number_controls, number_measurements)
	%PSYS convert rolmipvar to psys system description of matlab
	%	Input;:
	%		this:					instance
	%		number_controls:		number of controls if the vertices are not of type ltisys
	%		number_measurements:	number of measurements if the vertices are not of type ltisys
	%	Output:
	%		sys:	parameter dependent system description of rolmipvar
	if nargin <= 1
		number_controls = -1;
	else
		if ~isscalar(number_controls) || ~isnumeric(number_controls) || number_controls <= 0
			error('ROLMIP:psys', 'Number of conrtols must be a numeric scalar.');
		end
	end
	if nargin <= 2
		number_measurements = -1;
	else
		if ~isscalar(number_measurements) || ~isnumeric(number_measurements) || number_measurements <= 0
			error('ROLMIP:psys', 'Number of measurements must be a numeric scalar.');
		end
	end
	vertices = this.vertices;
	if get_degree(this) > 1
		% TODO: higher degrees might also be possible because only vertices are important for the uncertain area and the polynomial structure could be reduced to linear by dropping the dependence
		warning('ROLMIP:psys:convert', 'Polynomial systems of order larger than 1 can not be converted to ''psys'', uncertainty description will be conservative.');
	end
	m = size(this, 1);
	n = size(this, 2);
	if number_controls > 0 && number_measurements > 0
		if m < number_measurements || n < number_controls
			error('ROLMIP:psys', 'Number of measurements and controls must not exceed system dimension.');
		end
		if m - number_measurements ~= n - number_controls
			error('ROLMIP:psys', 'State matrix must be square.');
		end
		number_states = n - number_controls;
	else
		number_states = -1;
	end
	systems =cell(vertices, 1);
	e_i = 1:vertices;
	isltisys = false(vertices, 1);
	szSystems = NaN(vertices, 2);
	parfor ii = 1:vertices
		s = evalpar(this, {e_i == ii});
		%if size(s, 1) ~= m || size(s, 2) ~= n
		%	error('ROLMIP:add', 'Dimension of matrices must be consistent.');
		%end
		if isnan(s(end, end))
			s(end, end) = -Inf;
		end
		isltisys(ii, 1) = ltisys.isltisys(s);
		if ~isltisys(ii, 1) && number_states > 0
			s = ltisys(s(1:number_states, 1:number_states), s(1:number_states, number_states + (1:number_controls)), s(number_states + (1:number_measurements), 1:number_states), s(number_states + (1:number_measurements), number_states + (1:number_controls)));
			%if size(s, 1) ~= m || size(s, 2) ~= n
			%	error('ROLMIP:add', 'Dimension of matrices must be consistent.');
			%end
			isltisys(ii, 1) = true;
		end
		systems{ii, 1} = s;
		szSystems(ii, :) = [size(s, 1), size(s, 2)];
	end
	if ~all(isltisys(:))
		error('ROLMIP:psys', 'Vertices must be of type ''ltisys'' to convert to ''psys''.');
	end
	if ~all(szSystems(:, 1) == szSystems(1, 1)) || ~all(szSystems(:, 2) == szSystems(1, 2))
		error('ROLMIP:psys', 'Vertex systems must have the same dimension.');
	end
	sys = psys(cat(2, systems{:}));
end