function [Alpha_score, Alpha_pos, Convergence_curve, X_error, U_duan] = EOLSO(SearchAgents_no, Max_iter, lb, ub, dim, U, I, ocv, R0, model_RC)  
% Initialize SO algorithm parameters  
vec_flag = [1, -1]; % Position update direction  
Threshold = 0.25; % Food threshold  
Thresold2 = 0.6; % Temperature threshold  
C1 = 0.5;  % Food factor constant  
C2 = 0.05; % Constant 2  
C3 = 2;    % Constant 3  

% Initialize the snake population using Chebyshev  
P = Chebyshev(SearchAgents_no, dim, ub, lb);  

% Calculate initial fitness  
for i = 1:SearchAgents_no  
    [fitness(i), X_error, U_duan] = fobj_RC_R0(P(i,:), U, I, ocv, R0, model_RC);  
end  

[GYbest, gbest] = min(fitness);  
Xfood = P(gbest,:);  

% Divide population into male and female snakes  
Nm = round(SearchAgents_no / 2);  
Nf = SearchAgents_no - Nm;  
Xm = P(1:Nm,:);  
Xf = P(Nm+1:SearchAgents_no,:);  
fitness_m = fitness(1:Nm);  
fitness_f = fitness(Nm+1:SearchAgents_no);  

[fitnessBest_m, gbest1] = min(fitness_m);  
Xbest_m = Xm(gbest1,:);  
[fitnessBest_f, gbest2] = min(fitness_f);  
Xbest_f = Xf(gbest2,:);  

Convergence_curve = zeros(1, Max_iter);  
DOBL_SO_best = zeros(Max_iter, dim);  

% Initialize Xnewm and Xnewf  
Xnewm = zeros(Nm, dim);  
Xnewf = zeros(Nf, dim);  

for t = 1:Max_iter  
    DOBL_SO_best(t,:) = Xfood;  
    
    % Elite opposite learning  
    if mod(t, round(Max_iter/350)) == 0  
        [Xm, fitness_m] = Elite_Opposite(Xm, dim, Nm, @(x) fobj_RC_R0(x, U, I, ocv, R0, model_RC));  
        [Xf, fitness_f] = Elite_Opposite(Xf, dim, Nf, @(x) fobj_RC_R0(x, U, I, ocv, R0, model_RC));  
    end  
    
    Temp = rand * exp(-((t) / Max_iter));  
    Q = C1 * exp(((t - Max_iter) / Max_iter));  
    if Q > 1  
        Q = 1;  
    end  
    
    % Global exploration phase  
    if Q < Threshold  
        for i = 1:Nm  
            for j = 1:dim  
                rand_leader_index = floor(Nm*rand()+1);  
                X_randm = Xm(rand_leader_index, :);  
                flag_index = floor(2*rand()+1);  
                Flag = vec_flag(flag_index);  
                Am = exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));  
                Xnewm(i,j) = X_randm(j) + Flag*C2*Am*((ub(j)-lb(j))*rand+lb(j));  
            end  
        end  
        for i = 1:Nf  
            for j = 1:dim  
                rand_leader_index = floor(Nf*rand()+1);  
                X_randf = Xf(rand_leader_index, :);  
                flag_index = floor(2*rand()+1);  
                Flag = vec_flag(flag_index);  
                Af = exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));  
                Xnewf(i,j) = X_randf(j) + Flag*C2*Af*((ub(j)-lb(j))*rand+lb(j));  
            end  
        end  
    else  
        % Local development phase  
        if Temp > Thresold2  
            for i = 1:Nm  
                flag_index = floor(2*rand()+1);  
                Flag = vec_flag(flag_index);  
                for j = 1:dim  
                    Xnewm(i,j) = Xfood(j) + C3*Flag*Temp*rand*(Xfood(j)-Xm(i,j));  
                end  
            end  
            for i = 1:Nf  
                flag_index = floor(2*rand()+1);  
                Flag = vec_flag(flag_index);  
                for j = 1:dim  
                    Xnewf(i,j) = Xfood(j) + Flag*C3*Temp*rand*(Xfood(j)-Xf(i,j));  
                end  
            end  
        else  
            if rand > 0.6 % Fight mode  
                for i = 1:Nm  
                    for j = 1:dim  
                        FM = exp(-(fitnessBest_f)/(fitness_m(i)+eps));  
                        Xnewm(i,j) = Xm(i,j) + C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));  
                    end  
                end  
                for i = 1:Nf  
                    for j = 1:dim  
                        FF = exp(-(fitnessBest_m)/(fitness_f(i)+eps));  
                        Xnewf(i,j) = Xf(i,j) + C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));  
                    end  
                end  
            else % Mating mode  
                for i = 1:Nm  
                    for j = 1:dim  
                        Mm = exp(-fitness_f(i)/(fitness_m(i)+eps));  
                        Xnewm(i,j) = Xm(i,j) + C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));  
                    end  
                end  
                for i = 1:Nf  
                    for j = 1:dim  
                        Mf = exp(-fitness_m(i)/(fitness_f(i)+eps));  
                        Xnewf(i,j) = Xf(i,j) + C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));  
                    end  
                end  
            end  
        end  
    end  
    
    % Update positions and apply bounds  
    for j = 1:Nm  
        Xnewm(j,:) = bound(Xnewm(j,:), lb, ub);  
        [y, ~, ~] = fobj_RC_R0(Xnewm(j,:), U, I, ocv, R0, model_RC);  
        if y < fitness_m(j)  
            fitness_m(j) = y;  
            Xm(j,:) = Xnewm(j,:);  
        end  
    end  
    
    for j = 1:Nf  
        Xnewf(j,:) = bound(Xnewf(j,:), lb, ub);  
        [y, ~, ~] = fobj_RC_R0(Xnewf(j,:), U, I, ocv, R0, model_RC);  
        if y < fitness_f(j)  
            fitness_f(j) = y;  
            Xf(j,:) = Xnewf(j,:);  
        end  
    end  
    
    % Update best positions  
    [Ybest1, gbest1] = min(fitness_m);  
    [Ybest2, gbest2] = min(fitness_f);  
    
    if Ybest1 < fitnessBest_m  
        Xbest_m = Xm(gbest1,:);  
        fitnessBest_m = Ybest1;  
    end  
    if Ybest2 < fitnessBest_f  
        Xbest_f = Xf(gbest2,:);  
        fitnessBest_f = Ybest2;  
    end  
    
    % Record convergence  
    Convergence_curve(t) = min(Ybest1, Ybest2);  
    
    
    % Update global best  
    if fitnessBest_m < fitnessBest_f  
        GYbest = fitnessBest_m;  
        Xfood = Xbest_m;  
    else  
        GYbest = fitnessBest_f;  
        Xfood = Xbest_f;  
    end  
    
    % Display iteration info  
    %disp(['EOLSO iteration ' num2str(t) ', best fitness: ' num2str(GYbest)]);  
end  

Alpha_pos = Xfood;  
Alpha_score = GYbest;  
[~, X_error, U_duan] = fobj_RC_R0(Alpha_pos, U, I, ocv, R0, model_RC);  
end

function [x] = bound(x, l, u)  
    x = max(x, l);  
    x = min(x, u);  
end  
% 辅助函数  
function [SalpPositions,SalpFitness]=Elite_Opposite(SalpPositions,dim,N,fobj)
    for i=1:dim
        a=max(SalpPositions(:,i));
        b=min(SalpPositions(:,i));
        range(:,i)=[a;b];  
    end
    st=rand();
    for j = 1:dim  
        for k = 1:N
            Positionnew(k,j)=st*(sum(range(:,j)))-SalpPositions(k,j);
            if  Positionnew(k,j)>range(1,j)|| Positionnew(k,j)<range(2,j)
                Positionnew(k,j)=rand()*sum(range(:,j));
            end
        end   
    end
    total_Positions=[SalpPositions;Positionnew];
    for n=1:2*N
        total_Fitness(n)=fobj(total_Positions(n,:));
    end
    [~,index]=sort(total_Fitness);
    SalpPositions=total_Positions(index(1:N),:);
    SalpFitness=total_Fitness(index(1:N));
end


