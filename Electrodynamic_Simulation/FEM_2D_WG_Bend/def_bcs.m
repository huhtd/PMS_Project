% Assign BCS
% BC

% Better way
mode=1;
switch mode
    case 1
        % TE Modes: Ez
        for i=1:Nn
            if (st_no(i)==1 )
              for j=1:Nn
                A(i,j)=0.;
              end
              A(i,i)=1.;
            end
        end
    case 2
        % TM Modes: Hz
        return;
        for i=1:Nn
            if (st_no(i)==1 )
                for j=1:Nn
                    A(i,j)=0.;
                end
                A(i,i)=1.;
                % i
                break;
            end
        end
end

return;

% Worse way
delta=0.000000001;
for i=1:Nn
   if(abs(y_no(i))<delta)
      for j=1:Nn
          A(i,j)=0.;
          % B(i,j)=0.;
      end
      A(i,i)=1.e12;
      % B(i,i)=1.;
   end
   if(abs(y_no(i)-0.05)<delta)
      for j=1:Nn
          A(i,j)=0.;
          % B(i,j)=0.;
      end
      A(i,i)=1.e12;
      % B(i,i)=1.;
   end
   if(abs(x_no(i))<delta)
      for j=1:Nn
          A(i,j)=0.;
          % B(i,j)=0.;
      end
      A(i,i)=1.e12;
      % B(i,i)=1.;
   end
   if(abs(x_no(i)-0.1)<delta)
      for j=1:Nn
          A(i,j)=0.;
          % B(i,j)=0.;
      end
      A(i,i)=1.e12;
      % B(i,i)=1.;
   end
end
