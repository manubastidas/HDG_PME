function [Sol_F,Sol_V,Sol_U] = SolucionHDG(element,coordinate,...
    interioredge,GammaEdges,u_prev,u_itprev,parm,Timedt,Jac,Alocal,...
    A,B,C,D,E,G,R,posF,posU,posV,posFlujo,BASE)

nElement   = size(element,1);
nEdgeInt   = size(interioredge,1);
nEdgeGamma = size(GammaEdges,1);
nEdgeTotal = nEdgeInt+nEdgeGamma;

% FULL MATRICES
M = cell(nElement,1);
X = cell(nElement,1);
Y = cell(nElement,1);
% Inicialization - Indices
ipos = []; jpos = []; ipos2 = [];
Hval = []; Kval = [];

for elemento = 1:nElement
    
    % Coord (x;y)
    coord = coordinate(element(elemento,:),:)';
%     Afin  = coord(:,[2,3])-coord(:,[1,1]);
    
    %% RHS - F
    u_it1    = full(u_itprev(posU{elemento,1},:));
    u_prev1  = u_prev(posU{elemento,1});
    Llin     = parm*max(u_it1).^(parm-1);
    
    FF{elemento,1} =  matrizF(Jac(elemento,1),Alocal,u_prev1,u_it1,Llin,BASE,parm);
    
    %% LOCAL SOLVERS
    zero_aux0     = sparse(size(A{elemento,1},1),size(D{elemento,1},1));
    % Local matrix
    M{elemento,1} = [A{elemento,1} -B{elemento,1} zero_aux0;...
        B{elemento,1}' C{elemento,1} D{elemento,1};...
        zero_aux0' (1/Timedt)*D{elemento,1}' diag(-Llin.*diag(D{elemento,1}))];
    % Right hand side
    zero_aux      = sparse(size(C{elemento,1},1),size(G{elemento,1},2));
    D1{elemento,1} = [-E{elemento,1};G{elemento,1};zero_aux];
    D2{elemento,1} = [E{elemento,1};G{elemento,1};zero_aux];
    % ENSEMBLE
    X{elemento,1} = D2{elemento,1}'*(M{elemento,1}\D1{elemento,1})-R{elemento,1};
    Y{elemento,1} = (D2{elemento,1}'/M{elemento,1})*FF{elemento,1};
    
    % Optimization of the ensamble
    ipos     = [ipos;repmat(posFlujo{elemento,1}(:),size(posFlujo{elemento,1}(:),1),1)];
    jpos_aux = repmat(posFlujo{elemento,1}(:),1,size(posFlujo{elemento,1}(:),1))';
    jpos  = [jpos;  jpos_aux(:)];
    ipos2 = [ipos2; posFlujo{elemento,1}(:)];
    Hval = [Hval; X{elemento,1}(:)];
    Kval = [Kval; Y{elemento,1}];
end
H = sparse(ipos,jpos,Hval,BASE.cantFunciones2D*nEdgeTotal,BASE.cantFunciones2D*nEdgeTotal);
K = accumarray(ipos2,Kval);

%% Solution
SolFlujo = H\(-K);

%% Local solvers
% Initialization
posiF = []; posiV = []; posiU = [];
solF_val = []; solV_val = []; solU_val=[];

% Solution of the local solvers
parfor elemento = 1:nElement
    posflujoLocal  = posFlujo{elemento,1}(:);
    RHSlocal = D1{elemento,1}*SolFlujo(posflujoLocal,1)+FF{elemento,1};
    solLocal = M{elemento,1}\RHSlocal;
    
    posiF = [posiF;posF{elemento,1}'];
    posiV = [posiV;posV{elemento,1}'];
    posiU = [posiU;posU{elemento,1}'];
    
    solF_val = [solF_val;solLocal(1:2*BASE.cantFunciones)];
    solV_val = [solV_val;solLocal(2*BASE.cantFunciones+1:3*BASE.cantFunciones)];
    solU_val = [solU_val;solLocal(3*BASE.cantFunciones+1:end)];
end

% Rename the output
Sol_F = zeros(2*BASE.cantFunciones*nElement,1);
Sol_V = zeros(BASE.cantFunciones*nElement,1);
Sol_U = zeros(BASE.cantFunciones*nElement,1);

Sol_F(posiF,1)= solF_val;
Sol_V(posiV,1)= solV_val;
Sol_U(posiU,1)= solU_val;
