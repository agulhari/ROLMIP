function Z = mldivide(Y, X)
%MLDIVIDE (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9

  if(isa(X, 'rolmipvar'))
      if(isa(Y, 'rolmipvar'))
          error('ROLMIP:mrdivide', 'Division is not defined for variables of type ''rolmipvar''.');
      else
          Z = X;
          for i = 1:length(X.data)
            Z.data(i).value = Y\X.data(i).value;
            if (isscalar(Y))
                if (isa(Y,'sdpvar'))
                    Z.opcode{i} = ['<>', '\', X.opcode{i}];
                else
                    if (isreal(Y))
                        Z.opcode{i} = [num2str(Y), '#K1', '\', X.opcode{i}];
                    else %Imaginary
                        Z.opcode{i} = ['(', num2str(Y),')', '#K1', '\', X.opcode{i}];
                    end
                end
            else
                Z.opcode{i} = ['<>', '\', X.opcode{i}];
            end
          end
          if (isscalar(Y))
              if (isa(Y,'sdpvar'))
                  Z.label = ['<>', '\', X.label];
              else
                  if (isreal(Y))
                      Z.label = [num2str(Y), '\', X.label];
                  else %Imaginary
                      Z.label = ['(', num2str(Y), ')\', X.label];
                  end
              end
          else
              Z.label = ['<>', '\', X.label];
          end
      end
  elseif(isa(Y, 'rolmipvar'))
      error('ROLMIP:mrdivide', 'Division is not defined for variables of type ''rolmipvar''.');
  else
      % Why are you here?
      Z = Y\X;
  end
  % disp(['mtimes.m: ' Z.label ' has ' num2str(Z.vertices) ' vertices'])
  
end
