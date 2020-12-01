

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T: number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n: size of state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M: number of scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iter_stodcup: number of stodcup iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nb_iter _max: maximal number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xis{1,t}(k,j): scenario k component j for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Us{1,t}(k): value of U_t for scenario k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psis{1,t}(k): value of Psi_t for scenario k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tol: tolerance for the stopping criterion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upper_bounds _sddp(k) is the upper bound for SDDP iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower_bounds _sddp(k) is the lower bound for SDDP iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_sddp(k) is the CPU time needed to solve iteration k of the problem
% with SDDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upper_bounds _stodcup(k) is the upper bound for StoDCUP iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower_bounds _stodcup(k) is the lower bound for StoDCUP iteration k.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_stodcup(k) is the CPU time needed to solve iteration k of the problem
% with StoDCUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lower_bounds_sddp,upper_bounds_sddp,time_sddp]=inexact_sddp_quadratic(T,n,M,iter_inexact,nb_iter_max,xis,Us,Psis,tol,x0,probabilities,accuracies)

lower_bounds_sddp=[];
upper_bounds_sddp=[];
time_sddp=[];

subia=cell(1,T-1);
subja=cell(1,T-1);
valija=cell(1,T-1);

thetas=cell(1,T-1);

for t=1:(T-1)
    subia{1,t}=[subia{1,t},ones(1,n+1)];
    subja{1,t}=[subja{1,t},[2:1:n+2]];
    valija{1,t}=[valija{1,t},1,zeros(1,n)];
    thetas{1,t}= -(10^(10));
end

Cum_Probas=cell(1,T-1); 
for t=1:T-1
    Cum_Probas{1,t}=[0,cumsum(probabilities{1,t})];
end

End_Algo=1;
iter=1;
Costs=[];

maxe1=[];
maxe2=[];

while End_Algo
    iter
    %Forward pass
    total_cost=0;
    taux1=0;
    trial_states=cell(1,T-1);
    t=1;
    contloop=1;
    while ((t<=T)&&(contloop)) 
        if (t==1)
            Index=1;
        else
            Alea_Uniform=rand;
            [~,Index] = histc(Alea_Uniform,Cum_Probas{1,t-1});
            if (Alea_Uniform==1)
                Index=M;
            end
        end     
        tic
        xi=xis{1,t}(Index,:)';
        ur=Us{1,t}(Index);
        psir=Psis{1,t}(Index);
        clear prob;
        prob.qcsubk=[];
        prob.qcsubi=[];
        prob.qcsubj=[];
        prob.qcval=[];
        subi=[];
        subj=[];
        valij=[];  
        prob.c=[1,1,zeros(1,n)];
        if (t==T)
            prob.blx=[-inf;0;-100*ones(n,1)];
            prob.bux=[inf;0;100*ones(n,1)];
        else
            prob.blx=[-inf;-10^(10);-100*ones(n,1)];
            prob.bux=[inf;inf;100*ones(n,1)];
        end
        if (t==1)
            prob.buc=[-1-x0'*xi*xi'*x0;-ur;-4*n+psir;-1+psir;inf*ones(iter,1)];
        else
            if (t<T)
                prob.buc=[-1-trial_states{1,t-1}'*xi*xi'*trial_states{1,t-1};-ur;-4*n+psir;-1+psir;inf*ones(iter,1)];
            else
                prob.buc=[-1-trial_states{1,t-1}'*xi*xi'*trial_states{1,t-1};-ur;-4*n+psir;-1+psir];
            end
        end
        if (t<T)
            prob.blc=[-inf;-inf;-inf;-inf;thetas{1,t}];
        else
            prob.blc=[-inf;-inf;-inf;-inf];
        end
        %Constraint 1
        subi=[subi,ones(1,n+1)];
        subj=[subj,1,[3:1:n+2]];
        if (t==1)
            valij=[valij,-1,xi'-2*x0'*xi*xi'];
        else
            valij=[valij,-1,xi'-2*trial_states{1,t-1}'*xi*xi'];
        end
        for i=1:n
            for j=1:i
                prob.qcsubk=[prob.qcsubk,1];
                prob.qcsubi=[prob.qcsubi,2+i];
                prob.qcsubj=[prob.qcsubj,2+j];
                prob.qcval=[prob.qcval,2*xi(i)*xi(j)];
            end
        end
        
        %Constraint 2
        subi=[subi,2*ones(1,n+1)];
        subj=[subj,1,[3:1:n+2]];
        valij=[valij,-1,ones(1,n)];
        for i=1:n
            for j=1:i
                prob.qcsubk=[prob.qcsubk,2];
                prob.qcsubi=[prob.qcsubi,2+i];
                prob.qcsubj=[prob.qcsubj,2+j];
                prob.qcval=[prob.qcval,2*xi(i)*xi(j)];
            end
        end
        
        %Constraint 3
        subi=[subi,3*ones(1,n)];
        subj=[subj,[3:1:n+2]];
        valij=[valij,-8*ones(1,n)];
        for i=1:n
            prob.qcsubk=[prob.qcsubk,3];
            prob.qcsubi=[prob.qcsubi,2+i];
            prob.qcsubj=[prob.qcsubj,2+i];
            prob.qcval=[prob.qcval,8];
        end
        
        %Constraint 4
        subi=[subi,4*ones(1,n)];
        subj=[subj,[3:1:n+2]];
        valij=[valij,xi'];
        for i=1:n
            for j=1:i
                prob.qcsubk=[prob.qcsubk,4];
                prob.qcsubi=[prob.qcsubi,2+i];
                prob.qcsubj=[prob.qcsubj,2+j];
                prob.qcval=[prob.qcval,2*xi(i)*xi(j)];
            end
        end
        
        if (t<T)
            subi=[subi,subia{1,t}+4];
            subj=[subj,subja{1,t}];
            valij=[valij,valija{1,t}];
        end
        
        if (t==T)
            prob.a=sparse(subi,subj,valij,4,n+2);
        else
            prob.a=sparse(subi,subj,valij,4+iter,n+2);
        end
        if ((iter<=iter_inexact)&&(t>1))
            param.MSK_DPAR_INTPNT_TOL_REL_GAP=accuracies(iter);
            [~,res]=mosekopt('minimize echo(0)',prob,param);
        else
            [~,res]=mosekopt('minimize echo(0)',prob);
        end
        sol=res.sol.itr.xx;
        solsta=strcat('MSK_SOL_STA_', res.sol.itr.solsta);
        if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
            disp('Unfeasible primal problem F');
            contloop=0;
        elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
            disp('Primal infinite optimal value F');
            contloop=0;
        elseif (strcmp(solsta,'MSK_PRO_STA_UNKNOWN')==1)
            disp('Unknown status F');
            contloop=0;       
        else
            if (t<T)
                trial_states{1,t}=sol(3:n+2);
            end
        end
        if (contloop)
            if (t==1)
                zinf=sol(1)+sol(2);
            end
            total_cost=total_cost+sol(1);
        end
        taux1=taux1+toc;
        t=t+1;
    end
    
    if (contloop)
        Costs=[Costs;total_cost];
    end
    if (iter>=200)
        zsup=mean(Costs(iter-199:iter))+1.96*sqrt(var(Costs(iter-199:iter)))/sqrt(200);
        upper_bounds_sddp=[upper_bounds_sddp;zsup];
        lower_bounds_sddp=[lower_bounds_sddp;zinf];
    end
    % Backward pass
    t=T;
    slopes=zeros(T-1,n); 
    intercepts=zeros(T-1,1);
    
    while ((t>=2)&&(contloop))
        intercept=0;
        slope=zeros(n,1);
        j=1;
        while((j<=M)&&(contloop))
            tic
            xi=xis{1,t}(j,:)';
            ur=Us{1,t}(j);
            psir=Psis{1,t}(j);
            
            clear prob;
            prob.qcsubk=[];
            prob.qcsubi=[];
            prob.qcsubj=[];
            prob.qcval=[];
            subi=[];
            subj=[];
            valij=[];
            prob.c=[1,1,zeros(1,n)];
            if (t==T)
                prob.blx=[-inf;0;-100*ones(n,1)];
                prob.bux=[inf;0;100*ones(n,1)];
            else
                prob.blx=[-inf;-10^(10);-100*ones(n,1)];
                prob.bux=[inf;inf;100*ones(n,1)];
            end
            if (t==1)
                prob.buc=[-1-x0'*xi*xi'*x0;-ur;-4*n+psir;-1+psir;inf*ones(iter+1,1)];
            else
                if (t<T)
                    prob.buc=[-1-trial_states{1,t-1}'*xi*xi'*trial_states{1,t-1};-ur;-4*n+psir;-1+psir;inf*ones(iter+1,1)];
                else
                    prob.buc=[-1-trial_states{1,t-1}'*xi*xi'*trial_states{1,t-1};-ur;-4*n+psir;-1+psir];
                end
            end
            if (t<T)
                prob.blc=[-inf;-inf;-inf;-inf;thetas{1,t}];
            else
                prob.blc=[-inf;-inf;-inf;-inf];
            end
            
            %Constraint 1
            subi=[subi,ones(1,n+1)];
            subj=[subj,1,[3:1:n+2]];
            if (t==1)
                valij=[valij,-1,xi'-2*x0'*xi*xi'];
            else
                valij=[valij,-1,xi'-2*trial_states{1,t-1}'*xi*xi'];
            end
            for i=1:n
                for k=1:i
                    prob.qcsubk=[prob.qcsubk,1];
                    prob.qcsubi=[prob.qcsubi,2+i];
                    prob.qcsubj=[prob.qcsubj,2+k];
                    prob.qcval=[prob.qcval,2*xi(i)*xi(k)];
                end
            end
            
            %Constraint 2
            subi=[subi,2*ones(1,n+1)];
            subj=[subj,1,[3:1:n+2]];
            valij=[valij,-1,ones(1,n)];
            for i=1:n
                for k=1:i
                    prob.qcsubk=[prob.qcsubk,2];
                    prob.qcsubi=[prob.qcsubi,2+i];
                    prob.qcsubj=[prob.qcsubj,2+k];
                    prob.qcval=[prob.qcval,2*xi(i)*xi(k)];
                end
            end
            
            %Constraint 3
            subi=[subi,3*ones(1,n)];
            subj=[subj,[3:1:n+2]];
            valij=[valij,-8*ones(1,n)];
            for i=1:n
                prob.qcsubk=[prob.qcsubk,3];
                prob.qcsubi=[prob.qcsubi,2+i];
                prob.qcsubj=[prob.qcsubj,2+i];
                prob.qcval=[prob.qcval,8];
            end
            
            %Constraint 4
            subi=[subi,4*ones(1,n)];
            subj=[subj,[3:1:n+2]];
            valij=[valij,xi'];
            for i=1:n
                for k=1:i
                    prob.qcsubk=[prob.qcsubk,4];
                    prob.qcsubi=[prob.qcsubi,2+i];
                    prob.qcsubj=[prob.qcsubj,2+k];
                    prob.qcval=[prob.qcval,2*xi(i)*xi(k)];
                end
            end
            
            if (t<T)
                subi=[subi,subia{1,t}+4];
                subj=[subj,subja{1,t}];
                valij=[valij,valija{1,t}];
            end
            
            if (t==T)
                prob.a=sparse(subi,subj,valij,4,n+2);
            else
                prob.a=sparse(subi,subj,valij,5+iter,n+2);
            end
            if (iter<=iter_inexact)
                param.MSK_DPAR_INTPNT_TOL_REL_GAP=accuracies(iter);
                [~,res]=mosekopt('minimize echo(0)',prob,param);
            else
                [~,res]=mosekopt('minimize echo(0)',prob);
            end
            t
            iter
            if (iter<=iter_inexact)
                sol=res.sol.itr.xx;
                solsta=strcat('MSK_SOL_STA_', res.sol.itr.solsta);
                if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                    disp('Unfeasible primal problem B');
                    contloop=0;
                elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                    disp('Primal infinite optimal value B');
                    contloop=0;
                elseif (strcmp(solsta,'MSK_PRO_STA_UNKNOWN')==1)
                    disp('Unknown status B');
                    contloop=0;
                else
                    solsub=sol;
                    dual=res.sol.itr.suc;
                    psol=solsub(3:n+2);
                    auxc=zeros(1,n+2);
                    auxc(1)=1-dual(1);
                    auxc(2)=1;
                    auxc(3:n+2)=dual(1)*(2*xi*xi'*(psol-trial_states{1,t-1})+xi)';
                    lag=sol(1)+dual(1)*(psol'*xi*xi'*psol+psol'*(xi-2*xi*xi'*trial_states{1,t-1})-sol(1)+1+trial_states{1,t-1}'*xi*xi'*trial_states{1,t-1})-auxc*solsub;
                    slopeaux=dual(1)*2*xi*xi'*(trial_states{1,t-1}-psol);
                    slope=slope+slopeaux; 
                    intercept=intercept+lag-slopeaux'*trial_states{1,t-1};
                    
                    clear prob;
                    prob.qcsubk=[];
                    prob.qcsubi=[];
                    prob.qcsubj=[];
                    prob.qcval=[];
                    subi=[];
                    subj=[];
                    valij=[];
                    rac=(xi*xi')^(0.5);
                    if (t==T)
                        prob.blx=[-inf;0;-100*ones(n,1);-inf*ones(n,1)];
                        prob.bux=[inf;0;100*ones(n,1);inf*ones(n,1)];
                    else
                        prob.blx=[-inf;-10^(10);-100*ones(n,1);-inf*ones(n,1)];
                        prob.bux=[inf;inf;100*ones(n,1);inf*ones(n,1)];
                    end
                    if (t<T)
                        prob.buc=[-ur;-4*n+psir;-1+psir;inf*ones(iter+1,1);zeros(n,1)];
                    else
                        prob.buc=[-ur;-4*n+psir;-1+psir;zeros(n,1)];
                    end
                    if (t<T)
                        prob.blc=[-inf;-inf;-inf;thetas{1,t};zeros(n,1)];
                    else
                        prob.blc=[-inf;-inf;-inf;zeros(n,1)];
                    end
                    
                    %Constraint 2
                    subi=[subi,ones(1,n+1)];
                    subj=[subj,1,[3:1:n+2]];
                    valij=[valij,-1,ones(1,n)];
                    for i=1:n
                            prob.qcsubk=[prob.qcsubk,1];
                            prob.qcsubi=[prob.qcsubi,2+n+i];
                            prob.qcsubj=[prob.qcsubj,2+n+i];
                            prob.qcval=[prob.qcval,1];
                    end
                    
                    %Constraint 3
                    subi=[subi,2*ones(1,n)];
                    subj=[subj,[3:1:n+2]];
                    valij=[valij,-8*ones(1,n)];
                    for i=1:n
                        prob.qcsubk=[prob.qcsubk,2];
                        prob.qcsubi=[prob.qcsubi,2+i];
                        prob.qcsubj=[prob.qcsubj,2+i];
                        prob.qcval=[prob.qcval,8];
                    end
                    
                    %Constraint 4
                    subi=[subi,3*ones(1,n)];
                    subj=[subj,[3:1:n+2]];
                    valij=[valij,xi'];
                    for i=1:n
                            prob.qcsubk=[prob.qcsubk,3];
                            prob.qcsubi=[prob.qcsubi,2+n+i];
                            prob.qcsubj=[prob.qcsubj,2+n+i];
                            prob.qcval=[prob.qcval,1];
                    end
                    
                    if (t<T)
                        subi=[subi,subia{1,t}+3];
                        subj=[subj,subja{1,t}];
                        valij=[valij,valija{1,t}];
                        for i=1:n
                            subi=[subi,(4+iter+i)*ones(1,n+1)];
                            subj=[subj,2+n+i,[3:1:n+2]];
                            valij=[valij,1,-sqrt(2)*rac(i,:)];
                        end
                    else
                        for i=1:n
                            subi=[subi,(3+i)*ones(1,n+1)];
                            subj=[subj,2+n+i,[3:1:n+2]];
                            valij=[valij,1,-sqrt(2)*rac(i,:)];
                        end
                    end
                                        
                    if (t==T)
                        prob.a=sparse(subi,subj,valij,3+n,2*n+2);
                        nbc=3;
                    else
                        prob.a=sparse(subi,subj,valij,4+iter+n,2*n+2);
                        nbc=4+iter;
                    end
                    prob.c=[auxc,zeros(1,n)];
                    [~,res]=mosekopt('minimize echo(0)',prob);
                    solsta=strcat('MSK_SOL_STA_', res.sol.itr.solsta);
                    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                        disp('Unfeasible primal problem N');
                        contloop=0;
                    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                        disp('Primal infinite optimal value N');
                        contloop=0;
                    elseif (strcmp(solsta,'MSK_PRO_STA_UNKNOWN')==1)
                        disp('Unknown status N');
                        contloop=0;
                    else
                        solnoise=res.sol.itr.xx;
                        noiseterm=min([auxc*solnoise(1:n+2),auxc*solsub]);
                        noiseterm=noiseterm-auxc(1)*solsub(1)-auxc(3:n+2)*solsub(3:n+2);
                        intercept=intercept+noiseterm;
                    end
                end
                taux1=taux1+toc;
            else
                sol=res.sol.itr.xx;
                solsta=strcat('MSK_SOL_STA_', res.sol.itr.solsta);
                if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                    disp('Unfeasible primal problem');
                    contloop=0;
                elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                    disp('Primal infinite optimal value');
                    contloop=0;
                elseif (strcmp(solsta,'MSK_PRO_STA_UNKNOWN')==1)
                    disp('Unknown status');
                    contloop=0;
                else
                    dual=res.sol.itr.suc(1);
                    slopeaux=dual*2*xi*xi'*(trial_states{1,t-1}-sol(3:n+2));
                    slope=slope+slopeaux;
                    intercept=intercept+sol(1)+sol(2)-slopeaux'*trial_states{1,t-1};
                end
                taux1=taux1+toc;
            end
            j=j+1;
        end 
        if (contloop)
            subia{1,t-1}=[subia{1,t-1},(iter+1)*ones(1,n+1)];
            subja{1,t-1}=[subja{1,t-1},[2:1:n+2]];
            valija{1,t-1}=[valija{1,t-1},1,-slope'/M];
            thetas{1,t-1}=[thetas{1,t-1};intercept/M];
        else
            indexbreak=t;
        end
        t=t-1;
    end
    if (contloop==0)
        for t=(indexbreak+1):T
            subia{1,t-1}=subia{1,t-1}(1:iter*(n+1));
            subja{1,t-1}=subja{1,t-1}(1:iter*(n+1));
            valija{1,t-1}=valija{1,t-1}(1:iter*(n+1));
            thetas{1,t-1}=thetas{1,t-1}(1:iter);
        end
    end
    time_sddp=[time_sddp;taux1];
    if (contloop)
    End_Algo=abs((zsup-zinf)/zsup)>tol;
    if (iter>=nb_iter_max)
        End_Algo=0;
    end
    iter=iter+1;
    end
end





