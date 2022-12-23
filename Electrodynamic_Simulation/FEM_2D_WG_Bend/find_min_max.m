% Find min and max
max_field=abs(Vc(1));
for i=1:Nn
    if(abs(Vc(i))>max_field)
        max_field=abs(Vc(i));
    end
    %if(V(i)<min_field)
    %    min_field=V(i);
    %end
end
min_field=-max_field;
