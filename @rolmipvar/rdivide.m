function Z = rdivide(X, Y)
%RDIVIDE (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9

  if(isa(X, 'rolmipvar'))
      if(isa(Y, 'rolmipvar'))
          error('ROLMIP:rdivide', 'Division is not defined for variables of type ''rolmipvar''.');
      else
          Z = X;
          for i = 1:length(X.data)
            Z.data(i).value = X.data(i).value ./ Y;
            if (isscalar(Y))
                if (isa(Y,'sdpvar'))
                    Z.opcode{i} = [X.opcode{i}, './', '<>'];
                else
                    if (isreal(Y))
                        Z.opcode{i} = [X.opcode{i}, './', num2str(Y), '#K1'];
                    else %Imaginary
                        Z.opcode{i} = [X.opcode{i}, './(', num2str(Y),')', '#K1'];
                    end
                end
            else
                Z.opcode{i} = [X.opcode{i}, './', '<>'];
            end
          end
          if (isscalar(Y))
              if (isa(Y,'sdpvar'))
                  Z.label = [X.label, './', '<>'];
              else
                  if (isreal(Y))
                      Z.label = [X.label, './', num2str(Y)];
                  else %Imaginary
                      Z.label = [X.label, './(', num2str(Y),')'];
                  end
              end
          else
              Z.label = [X.label, './', '<>'];
          end
      end
  elseif(isa(Y, 'rolmipvar'))
      error('ROLMIP:mrdivide', 'Division is not defined for variables of type ''rolmipvar''.');
  else
      % Why are you here?
      Z = X./Y;
  end
  % disp(['mtimes.m: ' Z.label ' has ' num2str(Z.vertices) ' vertices'])
  
end
