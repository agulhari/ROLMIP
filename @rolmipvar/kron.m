function Z = kron(X, Y)
%KRON (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9

  if(isa(X, 'rolmipvar'))
      if(isa(Y, 'rolmipvar'))
          error('ROLMIP:kron', 'Kronecker product is not defined for variables of type ''rolmipvar''.');
      else
          Z = X;
          for i = 1:length(X.data)
            Z.data(i).value = kron(X.data(i).value, Y);
            if (isscalar(Y))
                if (isa(Y,'sdpvar'))
                    Z.opcode{i} = [X.opcode{i}, 'O', '<>'];
                else
                    if (isreal(Y))
                        Z.opcode{i} = [X.opcode{i}, 'O', num2str(Y), '#K1'];
                    else %Imaginary
                        Z.opcode{i} = [X.opcode{i}, 'O(', num2str(Y),')', '#K1'];
                    end
                end
            else
                Z.opcode{i} = [X.opcode{i}, 'O', '<>'];
            end
          end
          if (isscalar(Y))
              if (isa(Y,'sdpvar'))
                  Z.label = [X.label, 'O', '<>'];
              else
                  if (isreal(Y))
                      Z.label = [X.label, 'O', num2str(Y)];
                  else %Imaginary
                      Z.label = [X.label, 'O(', num2str(Y),')'];
                  end
              end
          else
              Z.label = [X.label, 'O', '<>'];
          end
      end
  elseif(isa(Y, 'rolmipvar'))
      Z = Y;
      for i = 1:length(Y.data)
        Z.data(i).value = kron(X, Y.data(i).value);
        if (isscalar(X))
            if (isa(X,'sdpvar'))
                Z.opcode{i} = ['<>', 'O', Y.opcode{i}];
            else
                if (isreal(X))
                    Z.opcode{i} = [num2str(X), '#K1', 'O', Y.opcode{i}];
                else %Imaginary
                    Z.opcode{i} = ['(',num2str(X), ')#K1O', Y.opcode{i}];
                end
            end
        else
            Z.opcode{i} = ['<>', 'O', Y.opcode{i}];
        end
      end
      if (isscalar(X))
          if (isa(X,'sdpvar'))
              Z.label = ['<>', 'O', Y.label];
          else
              if (isreal(X))
                  Z.label = [num2str(X), 'O', Y.label];
              else %Imaginary
                  Z.label = ['(',num2str(X),')O', Y.label];
              end
          end
      else
          Z.label = ['<>', 'O', Y.label];
      end
  else
      % Why are you here?
      Z = kron(X, Y);
  end
  % disp(['kron.m: ' Z.label ' has ' num2str(Z.vertices) ' vertices'])
  
end