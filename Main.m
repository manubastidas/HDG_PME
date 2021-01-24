
function Main(Archivo,Orden,mesh,Timedt,TimetnSteps,parm)

% global tol Time
load(Archivo)
tol = 1E-6;

% Initialization of the vector solutions
Solucion_u = zeros(cantFunciones*nElement,TimetnSteps);
Solucion_F = zeros(2*cantFunciones*nElement,TimetnSteps);
Solucion_v = zeros(cantFunciones*nElement,TimetnSteps);

% Inital solution
[Solucion_u(:,1)] = initial_u(nElement,coordinate,element,ALocal,Afin,BASE,parm);

residual  = zeros(0,TimetnSteps);
totalIter  = zeros(TimetnSteps,1);

% TIME LOOP
for kk=2:TimetnSteps+1
    
    it = 1;
    residual(1,kk) = inf;
    
    u_prev     = Solucion_u(:,kk-1);
    u_itprev   = u_prev;
    
    % Non-linear iterations
    while residual(it,kk)>tol && it<100
        
        % Non-linear solver
        [Sol_F,Sol_V,Sol_U] = SolucionHDG(element,coordinate,...
            interioredge,GammaEdges,u_prev,u_itprev,parm,Timedt,...
            Jac,ALocal,A,B,C,D,E,G,R,posF,posU,posV,posFlujo,BASE);
        
        % ITERATIONS!
        it = it+1;
        residual(it,kk)  = norm(Sol_U- u_itprev);
        u_itprev   = Sol_U;
        
        fprintf('\n Iter %i - %.4e',it-1,residual(it,kk))
    end
    totalIter(kk) = it;
    
    Solucion_F(:,kk)   = Sol_F;
    Solucion_v(:,kk)   = Sol_V;
    Solucion_u(:,kk)   = Sol_U;
    
    fprintf('\n Time %i of %i \n',kk-1,TimetnSteps)
    
    % Save the solution time kk
    Archivo_save0 = sprintf('Sol_O%i_mesh%i_m%i_t_%i',Orden,mesh,parm,kk);
    save(Archivo_save0);
end

%% SAVE THE SOLUTION
Archivo_save = sprintf('Sol_O%i_mesh%i_m%i',Orden,mesh,parm);
save(Archivo_save);

