classdef AHCMOEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained multi-objective optimization framework
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population1 = Problem.Initialization();           
            totalcon = size(Population1(1,1).con,2);          
            Set_cons = [];                                    
            Population2 = Problem.Initialization();           
            Fitness1    = CalFitness(Population1.objs,Population1.cons,Set_cons,true);          
            Fitness2    = CalFitness(Population2.objs,Population2.cons,Set_cons,false);         
            CONS = [Population1.cons;Population2.cons];
            CONS(CONS<0) = 0;
            VAR0 = max(sum(CONS,2));
            if VAR0 == 0
               VAR0 = 1;
            end
            X=0;
            Gen = 2;                                          
            threshold1 = 1e-3;                                
            threshold2 = 1e-2;                                
            Objvalues(1) = sum(sum(Population2.objs,1));      
            Objvalues(2) = 1000;
            Tr = Problem.N;                                   
            tau = 10;                                         
            Obj_matrix = [];                                  
            Con_matrix = [];                                  
            TR = [];                                          
            TR(1) = 1;
            TR(2) = 1;
            G = 100;                                          
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                 result = is_update(Objvalues,Gen,Population2,Problem.N,threshold1,Problem.M);  
                 %--------------------------------------------------------%   
                 if (result == 1) || (Tr <= threshold2)                                       
                      if length(Set_cons) ~= totalcon                                        
                          Cs = ACS_HT(Obj_matrix,Con_matrix,Population2,Set_cons);           
                          Set_cons = [Set_cons, Cs];                                         
                      end
                      Obj_matrix = [];                                                       
                      Con_matrix = [];                                                        
                      [Offspring1, Offspring2] = ACM_S(Problem,Population1, Population2, Fitness1, Fitness2);            
                 else       
                      [Offspring1, Offspring2] = ACM_S(Problem,Population1, Population2, Fitness1, Fitness2);
                 end 
%                if    Tr <= 0.01
%                             Offspring3 = IGM_DR(Problem,Population1,Population2,tau);                
%                else
                      Offspring3 = [];            
%                end
                 %--------------------------------------------------------%    
                      cp  = (-log(VAR0)-6)/log(1-0.5);
                      VAR = VAR0*(1-X)^cp;
                      [Population1, Fitness1] = main_task_EnvironmentalSelection([Population1,Offspring1,Offspring2,Offspring3],Problem.N,Set_cons,true,VAR);
                      tr = sum(ismember(Offspring2.decs, Population1.decs, 'rows'))/Problem.N; 
                      TR=[TR, tr];                                                             
                      if  length(TR) <= G + 1                                                  
                          Tr = 1;
                      else
                          Tr = mean(TR(Gen:-1:Gen-G));                                         
                      end
                      if  length(Set_cons) ~= totalcon
                          [Population2, Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2,Offspring3],Problem.N,Set_cons,false);   
                      else
                          [Population2, Fitness2] = EnvironmentalSelection1([Population2,Offspring1,Offspring2,Offspring3],Problem.N,Set_cons,false);   
                      end
                      Gen = Gen + 1;                                                           
                      Objvalues(Gen) = sum(sum(abs(Population2.objs),1));                    
                      Obj_information = Population2.objs;                                    
                      Con_information = Population2.cons;                                    
                      Obj_matrix=[Obj_matrix; mean(Obj_information)];                        
                      Con_matrix=[Con_matrix; mean(Con_information)];                        
                      X = X + 1/(Problem.maxFE/Problem.N);
            end
        end
    end
end 