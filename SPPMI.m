function [ SPPMI ] = createSPPMIMtx(G , k)
%% Creating SPPMI Matrix using graph G
% Calculating Degrees for each node
  nodeDegrees = sum(G);  
  nodeDegrees2=sum(G,2); 
  W = sum(nodeDegrees); 
  SPPMI = G;
% use a loop to calculate Wij*W/(di*dj)
  [col,row,weights] = find(G);
  for i = 1:length(col)
          score = log(weights(i) * W / nodeDegrees2(col(i)) / nodeDegrees(row(i))) - log(k);

            SPPMI(col(i),row(i)) = score;
         
  end


end

