function Offspring = IGM_DR(Problem,Population1,Population2,tau)
%-------------------------------------------------------------------------%

max_point = max(Population1.decs);      
min_point = min(Population2.decs);      
N_P=size(Population1.decs, 1);          
M_individuals = [];                     
H_individuals = [];                     
indices = randperm(N_P, tau);           
P1=Population1.decs;                    
P2=Population2.decs;                    
M_individuals = [M_individuals; P1(indices, :)];                     
H_individuals = [H_individuals; P2(indices, :)];                     
M = [];                                 
H = [];                                 
for i=1: tau
    d_c = M_individuals(i,:) - min_point;                            
    distance1 = sqrt(sum((M_individuals(i,:) - min_point) .^ 2));    
    d_f = max_point - H_individuals(i,:);                            
    distance2 = sqrt(sum(max_point - (H_individuals(i,:)) .^ 2));    
    num_sample = N_P/(2*tau);                                                                            
    for j = 1: num_sample
        Off_individual_M = M_individuals(i,:) - j/(num_sample + 1) * d_c;
        Off_individual_H = H_individuals(i,:) + j/(num_sample + 1) * d_f;
        M = [M; Off_individual_M];                                   
        H = [H; Off_individual_H];                                   
    end
end
    Offspring = [M; H];                                              
    Offspring = Problem.Evaluation(Offspring);                       
end