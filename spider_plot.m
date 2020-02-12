function spider_plot(MDL,numst)
%subplot(2,2,3)
%axes(mdl)
mdl = MDL;
indx = min(max(10, numst+1),numel(MDL));
plot(1:indx, MDL(1:indx), 's-k','LineWidth', 3,'MarkerSize', 6)%make sure the range is smaller than 10
%plot(MDL, 's-k','LineWidth', 3,'MarkerSize', 6)
y1 = min(MDL(1:indx)); y2 = max(MDL(1:indx));
d = y2-y1;
y1 = y1-d*0.2; y2 = y2+d*0.1;
ylim([y1,y2])
xlim([0.5, indx+0.5])
set(gca, 'XGrid','on')
ylabel('MDL value')
xlabel('# of States')
hold on
plot([1,1]*numst, [y1,y2], 'r--','LineWidth',3)
%plot([1,1]*numst,'r--','LineWidth',3)
hold off
%This must be the part I care about the most
hold off