% group = ones(size(LineWise));
data = data1;
Total = table2array(data);
group = ones(size(Total));
% c = 1;
groupx =[];
for i = 1:size(group,2)
%     r = mod(i,2);
%     if r==0 
%        c = c+0.5;
%     else
%         c = i;
%     end
    groupx = [groupx;i*group(:,i)];
end
    
legends = data.Properties.VariableNames;
boxplot(Total,groupx,'Labels',legends);
h = gca;
h.FontSize = 15;
h.FontWeight = 'bold';
h.XTickLabelRotation = 30;
title('Line 1');
xlabel(' Treatment');
ylabel('Angle (degrees)');
hold on
xCenter = 1:size(Total,2); 
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:size(Total,2)
    plot(rand(size(Total(:,i)))*spread -(spread/2) + xCenter(i), Total(:,i),'k*','linewidth', 2)
end
hold off