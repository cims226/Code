function Cs = ACS_HT(Obj_matrix, Con_matrix, Population2, Set_cons)
%-------------------------------------------------------------------------%
CV = sum(max(0,Population2.cons),2);            
fr=length(find(CV<=0))/size(Con_matrix,1);      
Analyze_Con = [];                               
cv = Population2.cons;                          
[m, n]=size(cv);                                
for i = 1:n
     ifr1 = length(find(cv(:,i) > 0))/m;        
     a_cv1 = mean(sum(max(0,cv(:,i)),2));       
     Analyze_Con=[Analyze_Con; ifr1 a_cv1];     
     if ifr1 == 0
        Set_cons =[Set_cons, i];                
     end
end
Set_cons = unique(Set_cons);                    
if  fr ==1
    Cs=[1,1:size(Con_matrix,1)];                      
elseif  length(Set_cons) ==  n
    Cs=[];
else
    All_cons = [1:size(Con_matrix,2)];                
    Un_cons = setdiff(All_cons, Set_cons);            
    %---------------------------------------------------------------------%
    numColsA = size(Obj_matrix, 2);
    numColsB = size(Un_cons, 2);
    C=[];                                             
    for i = 1:numColsA
        obj_con = [];                                 
        for j = 1:numColsB
            colA = Obj_matrix(:, i);                  
            colB = Con_matrix(:, Un_cons(j));         
            [~, p] = corrcoef(colA, colB);            
            if min(p) < 0.05                                   
               obj_con = [obj_con; i, Un_cons(j)];    
            end
        end
        if  ~isempty(obj_con)
            CL = obj_con(:,2)';                       
        else
            CL=[];
        end
        C = [C, CL];                                  
    end
        Cs = unique(C);                               
    %----------------------------------------------------------------------%  
    if  isempty(Cs)  
        Un_Analyze_Con = Analyze_Con(Un_cons,:);    
        distanceMatrix = pdist2(Un_Analyze_Con, Un_Analyze_Con, 'euclidean');    
        similarityMatrix = -distanceMatrix.^2;                                   
        % Run affinity propagation
        [idx, netsim, dpsim, expref] = affinityPropagation(similarityMatrix); 
        numClusters = numel(unique(idx));           
        ecv=[];                                     
        Clus={};                                    
        for k = 1:numClusters
            clusterPoints = Un_Analyze_Con(idx == k, :);         
            av_f = mean(clusterPoints,1);                        
            ecv(k) = av_f(1,1) * av_f(1,2);                      
            Clus{k} = Un_cons(idx == k);                         
        end
           [maxValue, maxIndex] = max(ecv);                      
           Cs = Clus{maxIndex};                                  
    end
end
end