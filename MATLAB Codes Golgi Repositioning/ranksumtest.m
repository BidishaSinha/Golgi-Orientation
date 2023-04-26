data = shapeindexS1;
var = cell(1,size(data,2));
for i = 1:size(data,2)
    var(i) = data.Properties.VariableNames(i);
end
P = ones(1,nchoosek(size(data,2),2));
c =1;
for i = 1:size(data,2)-1
    for j = i+1:size(data,2)
        x = table2array(data(:,i));
        y = table2array(data(:,j));
        P(c) = ranksum(x,y);
        c = c+1;
    end
end

