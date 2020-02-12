function [G, Ij, mj, split_tree] = clustering_GCP(eff, groups)
%% clusetering all the segments one by one up to only one
len = length(groups(1,:));% number of segments
mj = groups(2,:)-groups(1,:)+1;% the number of data points in each segment
Ij = zeros(1,len);% the averaged intensity of each segment
G = struct([]);% the structure to store all the clustering history
for i = 1:len
    Ij(i) = mean(eff(groups(1,i):groups(2,i)));
    G(1).g(i).gg = i;
end
[G, split_tree] = sub_clustering(Ij, mj, G);% clustering
end

%% 
% for each cycle, this sub function cluster two segments into one
function [G, split_tree] = sub_clustering(Ij, mj, G)
n = numel(G(end).g);%这里主要可以看出来有多少个segments
if n <= 1
    split_tree = [];
    return
else
    M = ones(n, n)*(-inf);% the merit matrix,这是一个含有n*n 负无穷的表格
    for i = 1:n-1
        for j = i+1:n
            % calculating the merit for all the possible pairs
            M(i, j) = llr_merit(Ij(i), Ij(j), mj(i), mj(j));%通过计算给这种merit matrix 每一个表格一种计算
        end
    end
    temp = [];
    while numel(G(end).g) > 1;
        numel(G(end).g)
        [a, b] = size(M);
        [~, q] = max(M(:));
        j = ceil(q/a);% find the column（应该是因为A(:)会把数据变为一列）
        i = q-(j-1)*a;% find the row
        %%
        if numel(G(end).g) <= 30%% to record the split tree %感觉这个值帮助分组分大小
            n = numel(G(end).g);
            if isempty(temp)
                temp = ones(1, n)*(n+0.5);
                split_tree = zeros(4, 3*n-1);
                split_tree(1:2, 2*n:end) = n+0.5;
                split_tree(3, 2*n:end) = Ij;
                split_tree(4, 2*n:end) = mj;
            end
            split_tree(1, 2*n-2:2*n-1) = n;
            split_tree(2, 2*n-2:2*n-1) = temp([i,j]);
            split_tree(3, 2*n-2:2*n-1) = Ij([i,j]);
            split_tree(4, 2*n-2:2*n-1) = mj([i,j]);
            temp(i) = n;
            temp(j) = [];
        end
        G(end+1).g = G(end).g;% initialize G(n+1) 在G 里面创建了另一个叫g的文件
        G(end).g(i).gg = [G(end).g(i).gg, G(end).g(j).gg];% group segments i and j together%感觉就是对着最大值的位置取值了吧
        %i是row，j是column，i那个点就变成[A,B]了
        G(end).g(j) = [];% free another space，这样的话，位置B 就什么都没有了，数据会上移一点。
        Ij(i) = (Ij(i)*mj(i)+Ij(j)*mj(j))/(mj(i)+mj(j));% group Ij
        Ij(j) = [];
        mj(i) = mj(i)+mj(j);% group mj
        mj(j) = [];
        %% update the M
        M(:, i) = -inf; M(i, :) = -inf;%基本上都是在合并值和删除值
        M(:, j) = []; M(j, :) = [];
        for k = 1:numel(G(end).g)% pairing with the new cluster i
            if k < i
                M(k, i) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));%应该和merit matrix 有关
            elseif k > i
                M(i, k) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));
            end% if
        end% for k
    end% while
    split_tree(1:4,1) = [1; 2; Ij; mj];
end% if n<=1
end% function

%% the merit function for Gaussian model (but can be more general)
function llr = llr_merit(I1, I2, m1, m2)
I = (I1*m1 + I2*m2)/(m1+m2);
llr = (m1+m2)*I^2-m1*I1^2-m2*I2^2;
end