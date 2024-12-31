function [Offspring1, Offspring2] = ACM_S(Problem, Population1, Population2, Fitness1, Fitness2)
 %------------------------------------------------------------------------%   
     Offspring1={};
     Offspring2={};
     O1 = [];
     O2 = [];
     P1 = Population1.decs;
     P2 = Population2.decs;
     jsd_value = jensen_shannon_divergence(Population1.decs, Population2.decs);   
     if isreal(jsd_value)
        C_n = floor(Problem.N * jsd_value);     
     else
        C_n = 0;
     end
     max_point = max(Population1.decs);               
     min_point = min(Population2.decs);               
     %--------------------------------------------------------------------%
     for i=1:C_n                                       
         %----------------------------------------------------------------%
         distances1 = vecnorm(P2 - P1(i, :), 2, 2); 
         [~, idx] = sort(distances1); 
         refpoint = 1/3 * (min_point + P2(idx(1),:) + P2(idx(2),:));
         o1 = P1(1,:) + rand * (refpoint - P1(1,:));
         O1 = [O1; o1];
         %----------------------------------------------------------------%
         distances2 = vecnorm(P1 - P2(i, :), 2, 2); 
         [~, idx] = sort(distances2); 
         refpoint = 1/3 * (max_point + P1(idx(1),:) + P1(idx(2),:));
         o2 = P2(1,:) + rand * (refpoint - P2(1,:));
         O2 = [O2; o2];
     end
     if ~isempty(O1)  
        Offspring1 = Problem.Evaluation(O1);
     end
     if ~isempty(O2)
        Offspring2 = Problem.Evaluation(O2);
     end
     %--------------------------------------------------------------------%
     MatingPool1 = TournamentSelection(2,2*(Problem.N-C_n),Fitness1);
     MatingPool2 = TournamentSelection(2,2*(Problem.N-C_n),Fitness2);
     O3 = OperatorDE(Problem,Population1(C_n+1:Problem.N),Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
     O4 = OperatorDE(Problem,Population2(C_n+1:Problem.N),Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
     Offspring1 = [Offspring1, O3];
     Offspring2 = [Offspring2, O4];
end