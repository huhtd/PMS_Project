% Assemble the matrix

A=sparse(Nn,Nn);    % Sparse system matrix
% A=zeros(Nn,Nn);         % Dense system matrix
b=zeros(Nn,1);        % RHS
% B=sparse(Nn,Nn);
%for i=1:Nn
%    for j=1:Nn
%        A(i,j)=0.;
%        B(i,j)=0.;
%    end
%end

for ie=1:Ne
    % ke=int_cont_es(i);
    int_cont_wp;
    for k=1:3
        for l=1:3
            A(el_no(ie,k),el_no(ie,l))=A(el_no(ie,k),el_no(ie,l))+ae(k,l);
        end
    end
end
