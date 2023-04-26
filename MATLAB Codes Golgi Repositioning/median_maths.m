data = medians;
div = ones(1,nchoosek(size(data,2),2));
dif = div;
c = 1;
for i = 1:size(data,2)-1
    for j = i+1:size(data,2)
        div(c) = data(i)/data(j);
        dif(c) = abs(data(i)-data(j));
        c = c+1;
    end
end
