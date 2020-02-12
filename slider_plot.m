function [pl_r,pl_i]=slider_plot(eff,eff_fit,current_state)
%This must be the part I care about the most
try
    right = numel(eff);
catch ME
    return
end
if right <= 100
    return
end
%left = get(handles.leftpos,'Value');
%left = max(1,ceil(left*(right-100)));
%ratio = get(handles.range,'Value');
%rg = exp((log(right-left+1)-log(100))*ratio+log(100));
%rg = round(rg);
%rg=length(eff);
%right = min(left+rg,numel(handles.eff));

%axes(handles.traces)
x=1:length(eff(:));
numst = current_state;
pl_r=plot(x*0.1,eff(:), ':','LineWidth',1.25)
hold on
pl_i=plot(x*0.1,eff_fit(numst,:), 'r','LineWidth',2)
[eby, ebx] = hist(eff, 50);
y1 = min(ebx)-0.01; y2 = min(1.1,max(ebx)+0.01);
%ylim([y1,y2]);
%xlim([left,right]);
xlabel('Exp Time (s)')
ylabel('Counts (10 ms) ')
grid on
hold off