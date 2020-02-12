%The data source need to be vertical matrix,absolute time is in nanosecond
%Binsize is also in nanosecond
function timetrace=Gen_timetrace(datasource,binsize)

datasource_edge=[datasource(1,1):binsize:datasource(end,1)];
[N,edge]=histcounts(datasource,datasource_edge);
timetrace= cat(2,transpose(edge(1,2:end)),transpose(N));

%%%check point%%%
%scatter(edge(1,2:end),N);
%%%check point
    
end