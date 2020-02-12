function DandABS=GetDandABS(datasource,varargin);
for i=1:nargin-1
    if ischar(varargin{i})
        if strcmp(varargin{i},'M')
            try
            datasource=datasource(find(datasource(:,2)==3,1):end,:);
            catch
            error('This file do not have a marker.')    
            end
            stringvar=i;           
        end
    end
end
varargin(stringvar)=[];
switch nargin-2
    case 1
        photon=datasource(datasource(:,2)==9,:);
        DandABS=photon(photon(:,3)==varargin{1},5:6);
    case 2
        photon=datasource(datasource(:,2)==9,:);
        DandABS=photon(photon(:,3)==varargin{1} | photon(:,3)==varargin{2},5:6);
    otherwise
        error('More than 2 channel is not considered yet')
        
end



end