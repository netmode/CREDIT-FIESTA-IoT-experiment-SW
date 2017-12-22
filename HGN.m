function [modules] = HGN(adj,k,batchSize,o,l,x,t,r,e)
% Inputs:   
%           adj: the adjacency matrix of the graph
%           k: the required number of communities
%           batchSize: the size of the Batch (up to how many edges will be
%           removed before a new embedding)
%           The next parameters are used for the rigel embedding
%           o: path for the output file
%           l: path for the Dist file
%           x: the number of dimensions
%           t: input file path
%           r: path for the landmarks file
%           e: curvature of hyperbolic space
% Outputs:
%           modules: cell array 1 x k . Each cell contains the nodes of the
%           corresponding module.
% Author:   Konstantinos Tsitseklis           


n=size(adj,1);
modules{1}=[1:n];

curr_mod=1;

adj_temp=adj;

while length(modules)<k
    length(modules)
    u = size(adj_temp,1);
    L = 6;
    b = 2;
    i = 4;
    
    coordinatesMatrix = embed(adj_temp,o,l,num2str(b),num2str(x),t,r,num2str(u),num2str(e),num2str(L),num2str(i));

    iterations = size(adj_temp,1);
    w = hyperEBC(adj_temp,coordinatesMatrix,iterations);
    p = 1;
    while ((p<=2*batchSize) && (isconnected(adj_temp)))
        adj_temp(w(p,1),w(p,2)) = 0;
        adj_temp(w(p,2),w(p,1)) = 0;
        p = p+2;
    end

    if isconnected(adj_temp); continue; end % keep on removing edges

    comp_mat = find_conn_comp(adj_temp);
    
    for c=1:length(comp_mat); comp_mat{c};end
    
    for c=1:length(comp_mat); modules{length(modules)+1}=modules{curr_mod}(comp_mat{c}); end
    
    % remove "now" disconnected component (curr_mod) from modules
    modules{curr_mod}=modules{1};
    modules=modules(2:length(modules));
    
    modL=[];
    for j=1:length(modules); modL(j)=length(modules{j}); end
    [maxL,indL]=max(modL);
    curr_mod=indL;
    
    adj_temp=subgraph(adj,modules{indL});
end 

end