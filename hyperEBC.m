function [HEBC] = hyperEBC( adjMatrix, coordinatesMatrix, src )
% ALGORITHM FOR HYPERBOLIC BETWEENNESS CENTRALITY
% Inputs:
%           AdjMatrix: The adjacency matrix of the network
%           coordinatesMatrix: A n x m matrix in which every row represents
%           a node and every column the coordinate of the node in that
%           dimension (n nodes, m dimensions total) 
% Outputs: 
%           EBC: A |E|x 3 matrix. Each row represents an edge. Column 1
%           represents the starting and column 2 the finishing point (node)                     
%           of the edge. Column 3 is the EBC value of the edge specified 
%           by the previous columns.
%           Authors: Konstantinos Tsitseklis, Konstantinos Sotiropoulos.
%% Pre-Processing of Adjacency Matrix

%%
edgesNumber = numedges(adjMatrix);
dim = size(coordinatesMatrix,2);
% jm = java.util.HashMap;
tic
HEBC = zeros(edgesNumber*2,3);

% k=1;
% for i=1:size(adjMatrix,1)
%     for j=1:size(adjMatrix,1)
%         if (adjMatrix(i,j)==1)
%             HEBC(k,1) = i;
%             HEBC(k,2) = j;
%             s1 = num2str(i);
%             s2 = num2str(j);
%             s3 = strcat(s1,'-',s2);
%             jm.put(s3,0);
%             k = k+1;
%         end
%     end
% end


%%
nodesNumber = size( adjMatrix, 2 );
pSize = max(sum(adjMatrix));
testMatrix = zeros(nodesNumber, pSize);
indexMatrix = zeros(nodesNumber,1);
for i=1:nodesNumber
    tempArray = find(adjMatrix(i,:));
    for j=1:length( tempArray )
        testMatrix(i,j) = tempArray(j);
    end
    indexMatrix(i) = length( tempArray );
end

clear  tempArray

nodesNumber = ( size( coordinatesMatrix, 1 ) );

%%
temp = zeros( nodesNumber );
%%
for destination=1:src
    indexP = zeros( nodesNumber,1 );
    
    sigma = zeros( nodesNumber, 1);
    sigma( destination ) =1;
    distances = zeros(1,nodesNumber);
    %% STAGE 1 - TOPOLOGICAL SORT
    dst = coordinatesMatrix(destination,:);
    
    % a faster way to calculate nodes distance
    ysum = 1; %Sxi^2
    for j=2:dim
        ysum = ysum + dst(j).^2;
    end
    for vertex=1:nodesNumber
        %if ( vertex~=destination )
        xsum = 1; %Syi^2
        xysum = 0;
        
        for j=2:dim
            xsum = xsum + coordinatesMatrix(vertex,j).^2;
            xysum = xysum + coordinatesMatrix(vertex,j)*dst(j);
        end
        
        t = sqrt( ysum*xsum ) - xysum;
        dist = acosh(t);
        distances(vertex) = dist;
        %end
    end

    [~,DAG] = sort(distances, 'descend');
    %% PART 2 - number of greedy paths between nodes and destination
    for i=nodesNumber:-1:1
        v = DAG(i);
        for j=1:indexMatrix(v)
            w = testMatrix(v,j);
            if (( distances(w) > distances(v) + 0.3))
                sigma(w) = sigma(w)+sigma(v);
                indexP(w) = indexP(w)+1;
                P(w,indexP(w)) = v;
            end
        end
    end
    
    %% PART 3
    delta = zeros( nodesNumber, 1 );
    % S returns vertices in order of non-increasing distance from s
    for node=1:nodesNumber
        % pop w<-S
        w = DAG( node );
        if sigma(w)>0
            for j=1:indexP(w)
                v = P(w,j);
                %%
                c = ( sigma(v)/sigma(w) )*(1 + delta(w)) ;
                temp(w,v) = temp(w,v) + c;
                temp(v,w) = temp(v,w) + c;
%                 s1 = num2str(w);
%                 s2 = num2str(v);
%                 s3 = strcat(s1,'-',s2);
%                 hebc2 = jm.get(s3);
%                 hebc2 = hebc2 + c;
%                 jm.put(s3,hebc2);
%                 s3 = strcat(s2,'-',s1);
%                 jm.put(s3,hebc2);
                %%
                delta(v) = delta(v)+c ;
            end
        end
    end  
end
% HEBC = HEBC/edgesNumber;

toc
% sort edges according to their HEBC value in decreasing order
k=1;
for i=1:nodesNumber
    for j=1:nodesNumber
        if (adjMatrix(i,j)==1)
%             s1 = num2str(i);
%             s2 = num2str(j);
%             s3 = strcat(s1,'-',s2);
            HEBC(k,1) = i;
            HEBC(k,2) = j;
%             HEBC(k,3) = jm.get(s3);
            HEBC(k,3) = temp(i,j);
            k=k+1;
        end
    end
end

[~, E] = sort(HEBC(:,3),'descend');
HEBC = HEBC(E,:);

end