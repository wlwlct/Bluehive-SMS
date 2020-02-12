%%
function state_plot(numst,eff,eff_fit,MDL,current_state)
%This should be the point of findingding state
%axes(handles.histos)
[eby, ebx] = hist(eff, 50);
y1 = min(ebx)-0.01; y2 = min(1.1,max(ebx)+0.01);
barh(ebx, eby);
%ylim([y1,y2]);
%set(gca, 'Ytick',[],'Xtick',[]);
xlabel('Relative Population')
hold on
id = unique(eff_fit(numst,:));
if numel(id)==1
    hst = numel(eff);
else
    hst = hist(eff_fit(numst,:), id);
end
hst = hst./round(max(hst)/max(eby));
plot([zeros(size(hst));hst],[id;id], 'r','LineWidth',5)
plot(hst, id, 'r.', 'MarkerSize', 10)
hold off

%This should be th part of choosing which is the right data
%slider_plot(eff,eff_fit,current_state,MDL)
