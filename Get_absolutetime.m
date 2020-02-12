function [data1,data2]=Get_absolutetime(PTU_raw,channel1,varargin)
switch nargin
    case 4
        if strcmp(varargin{2},'M')
        try
        PTU_raw=PTU_raw(find(PTU_raw(:,2)==3,1):end,:);
        catch
        error('This file do not have a marker.')
        end
        end
        photon=PTU_raw(PTU_raw(:,2)==9,:);
        data1=photon(photon(:,3)==channel1,5);
        data2=photon(photon(:,3)==varargin{1},5);
    case 3 
        if strcmp(varargin,'M')
        try
        PTU_raw=PTU_raw(find(PTU_raw(:,2)==3,1):end,:);
        catch
        error('This file do not have a marker.')    
        end
        photon=PTU_raw(PTU_raw(:,2)==9,:);
        data1=photon(photon(:,3)==channel1,5);
        data2=[];
        else
        data1=photon(photon(:,3)==channel1,5);
        data2=photon(photon(:,3)==varargin,5);
        end
    otherwise
        error('Not correct input');
end
end