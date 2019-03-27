function Z = times(X, Y)
%TIMES (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9

  if(isa(X, 'rolmipvar'))
      if(isa(Y, 'rolmipvar'))
		  % TODO: implement it?
          error('ROLMIP:times', 'Element wise product of two rolmipvars is not defined.');
      else
          Z = X;
          for i = 1:length(X.data)
            Z.data(i).value = X.data(i).value .* Y;
            if (isscalar(Y))
                if (isa(Y,'sdpvar'))
                    Z.opcode{i} = [X.opcode{i}, '.*', '<>'];
                else
                    if (isreal(Y))
                        Z.opcode{i} = [X.opcode{i}, '.*', num2str(Y), '#K1'];
                    else %Imaginary
                        Z.opcode{i} = [X.opcode{i}, '.*(', num2str(Y),')', '#K1'];
                    end
                end
            else
                Z.opcode{i} = [X.opcode{i}, '.*', '<>'];
            end
          end
          if (isscalar(Y))
              if (isa(Y,'sdpvar'))
                  Z.label = [X.label, '.*', '<>'];
              else
                  if (isreal(Y))
                      Z.label = [X.label, '.*', num2str(Y)];
                  else %Imaginary
                      Z.label = [X.label, '.*(', num2str(Y),')'];
                  end
              end
          else
              Z.label = [X.label, '.*', '<>'];
          end
      end
  elseif(isa(Y, 'rolmipvar'))
      Z = Y;
      for i = 1:length(Y.data)
        Z.data(i).value = X .* Y.data(i).value;
        if (isscalar(X))
            if (isa(X,'sdpvar'))
                Z.opcode{i} = ['<>', '.*', Y.opcode{i}];
            else
                if (isreal(X))
                    Z.opcode{i} = [num2str(X), '#K1', '.*', Y.opcode{i}];
                else %Imaginary
                    Z.opcode{i} = ['(',num2str(X), ')#K1.*', Y.opcode{i}];
                end
            end
        else
            Z.opcode{i} = ['<>', '.*', Y.opcode{i}];
        end
      end
      if (isscalar(X))
          if (isa(X,'sdpvar'))
              Z.label = ['<>','.*', Y.label];
          else
              if (isreal(X))
                  Z.label = [num2str(X),'.*', Y.label];
              else %Imaginary
                  Z.label = ['(',num2str(X),').*', Y.label];
              end
          end
      else
          Z.label = ['<>','.*', Y.label];
      end
  else
      % Why are you here?
      Z = X.*Y;
  end
  % disp(['mtimes.m: ' Z.label ' has ' num2str(Z.vertices) ' vertices'])
  
end
