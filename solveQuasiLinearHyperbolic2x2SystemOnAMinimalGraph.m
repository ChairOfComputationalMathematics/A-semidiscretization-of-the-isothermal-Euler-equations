function [tspan, sol] = solveQuasiLinearHyperbolic2x2SystemOnAMinimalGraph()
%% Explanation
% Solving a quasi-linear system of transport equations with a source term
% on an easy graph
%              pL(t)->    p      p
%                     (0)--->(1)---(2)
%                         q      q    -> bR(t)
%
%           [ r11 ]   + [ lambda_11     0     ] [ r11 ]   = [ F11(t,x,R) ]
%           [ r12 ]_t   [    0      lambda_12 ] [ r12 ]_x = [ F12(t,x,R) ] 
%               
%           [ r21 ]   + [ lambda_21     0     ] [ r21 ]   = [ F21(t,x,R) ]
%           [ r22 ]_t + [    0      lambda_22 ] [ r22 ]_x = [ F22(t,x,R) ] 
%
% with boundary conditions (r_ijk(t) = r_ij(t,x_k))
%           (r111 + r121)(t) = 2log(pL(t))
%           (r21n+1 + r22n+1)(t) + 2log((r21n+1 - r22n+1)(t)) = 2log(abR(t))
% and coupling conditions:
%           (compressor property:)
%           (r211 + r221)(t) - (r11n+1 + r12n+1)(t) = 2log(u(t))
%           (conservation of mass:)
%           (r211 - r221)(t) - (r11n+1 - r12n+1)(t) = 0 WRONG!!! INSERT u*p!!!
% and initial condition
%
% on a domain [0,T] x [0,L] for both edges.
%

%% Initializations

T = 10;
L = 100;

a = 340;
lambdaFric = .15;
diam = .1;

% be careful, that C^0-compatibility condition are fulfilled for boundary
% and initial data
pressureLeft = @(t) 1e5;
dtPressureLeft = @(t) 0;
flowRight = @(t) 30;
dtFlowRight = @(t) 0;
control = @(t) 1.5;
dtcontrol = @(t) 0;
initialPressure1 = @(x) 1e5;
initialPressure2 = @(x) control(0)*1e5;
initialFlow1 = @(x) 30;
initialFlow2 = @(x) 30;
initialr11 = @(x) log(initialPressure1(x)) + a*initialFlow1(x)./initialPressure1(x);
initialr12 = @(x) log(initialPressure1(x)) - a*initialFlow1(x)./initialPressure1(x);
initialr21 = @(x) log(initialPressure2(x)) + a*initialFlow2(x)./initialPressure2(x);
initialr22 = @(x) log(initialPressure2(x)) - a*initialFlow2(x)./initialPressure2(x);

n = 100;
lspan = linspace(0,L,n+1);
deltax = L/n;

% mass matrix
MassMatrix = speye(4*n+4);
MassMatrix(2*n+2,n+1) = -1;
MassMatrix(2*n+2,2*n+2) = -1;
MassMatrix(2*n+2,2*n+3) = 1;
MassMatrix(2*n+2,3*n+4) = 1;
MassMatrix(2*n+3,n+1) = -1;
MassMatrix(2*n+3,2*n+2) = 1;

% complete right hand side dt R(t) = RHS(t,R(t))

RHS = @(t,R) inv(getMassMatrix(t,R)) * (getF(t,R) - getFDM(t,R)*R);

%% Crank-Nicolson method to solve remaining ODE

% initializations
numit = 1000;
tspan = linspace(0,T,numit);
deltat = T/(numit-1);
theta = .5;

sol = zeros(numit,4*n+4);

% initial data
y0 = zeros(1,4*n+4);
y0(1,1:n+1) = initialr11(lspan);
y0(1,n+2:2*n+2) = initialr12(lspan);
y0(1,2*n+3:3*n+3) = initialr21(lspan);
y0(1,3*n+4:4*n+4) = initialr22(lspan);
sol(1,:) = y0;

% Crank-Nicolson Method
for i = 2:numit         % loop over time steps
    
    t = tspan(i);
    yn = sol(i-1,:);    % last step solution
    yNew = yn;          % first estimate for current step solution
    
    while 1             % inner loop for better estimates
        
        yTheta = theta*yNew + (1-theta)*yn;     % mix between old and 
                                                % new solution
        
        %%% get MassMatrix M
        % boundary conditions
        %MassMatrix(1,1) = 1;
        MassMatrix(1,n+2) = 1;
        MassMatrix(end,3*n+3) = 1 + 2/(yTheta(3*n+3) - yTheta(4*n+4));
        MassMatrix(end,end) = 1 - 2/(yTheta(3*n+3) - yTheta(4*n+4));
        % coupling conditions
        %MassMatrix(2*n+2,n+1) = -1;
        %MassMatrix(2*n+2,2*n+2) = -1;
        %MassMatrix(2*n+2,2*n+3) = 1;
        %MassMatrix(2*n+2,3*n+4) = 1;
        %MassMatrix(2*n+3,n+1) = -1;
        %MassMatrix(2*n+3,2*n+2) = 1;
        MassMatrix(2*n+3,2*n+3) = control(t);
        MassMatrix(2*n+3,3*n+4) = -control(t);
        
        %%% get Lambda
        Lambda = .5*a*sparse(1:4*n+4,1:4*n+4,[yTheta(1:n+1)' - yTheta(n+2:2*n+2)' + 2; ...
                                              yTheta(1:n+1)' - yTheta(n+2:2*n+2)' - 2; ...
                                              yTheta(2*n+3:3*n+3)' - yTheta(3*n+4:4*n+4)' + 2; ...
                                              yTheta(2*n+3:3*n+3)' - yTheta(3*n+4:4*n+4)' - 2]);
        
        %%% get finite difference matrix
        % initialization
        V = ones(n, 1);
        D = diag(V, 1) - diag(V, -1);
        D(1,1:2) = [-2 2];
        D(n+1,n:n+1) = [-2 2];
        V = zeros(n+1);
        FDM = [D V V V; V D V V; V V D V; V V V D];
        FDM = sparse(FDM);
        FDM = Lambda * 2/deltax * FDM;
        % boundary conditions
        FDM(1,1:2) = [0 0];
        FDM(end,end-1:end) = [0 0];
        % coupling conditions
        FDM(2*n+2,2*n+1:2*n+2) = [0 0];
        FDM(2*n+3,2*n+3) = dtcontrol(t);
        FDM(2*n+3,2*n+4) = 0;
        FDM(2*n+3,3*n+4) = -dtcontrol(t);
        
        %%% get F
        % initialization
        F = zeros(4*n+4,1);
        F(1:n+1) = -1/8*a*lambdaFric/diam * (yTheta(1:n+1) - yTheta(n+2:2*n+2)).^2;
        F(n+2:2*n+2) = 1/8*a*lambdaFric/diam * (yTheta(1:n+1) - yTheta(n+2:2*n+2)).^2;
        F(2*n+3:3*n+3) = -1/8*a*lambdaFric/diam * (yTheta(2*n+3:3*n+3) - yTheta(3*n+4:4*n+4)).^2;
        F(3*n+4:4*n+4) = 1/8*a*lambdaFric/diam * (yTheta(2*n+3:3*n+3) - yTheta(3*n+4:4*n+4)).^2;
        % boundary conditions
        F(1) = 2/pressureLeft(t) * dtPressureLeft(t);
        F(end) = 2/flowRight(t) * dtFlowRight(t);
        % coupling conditions
        F(2*n+2) = 1/control(t) * dtcontrol(t);
        F(2*n+3) = 0;
        
        %yTmp = yn' + deltat*RHS(t,yTheta');     % better estimate for 
                                                % current step solution
        %yTmp = (getMassMatrix(t,yTheta) + deltat*theta*getFDM(t,yTheta))\((getMassMatrix(t,yTheta) - deltat*(1-theta)*getFDM(t,yTheta))*yn' + deltat*getF(t,yTheta));                                      
        yTmp = (MassMatrix + deltat*theta*FDM)\((MassMatrix - deltat*(1-theta)*FDM)*yn' + deltat*F);
        yTmp = yTmp';
        str = ['Error in step nr. ' num2str(i) ': ' num2str(norm(yTmp - yNew))];
        %disp(str);
        if norm(yTmp - yNew) < sqrt(sqrt(eps))        % exit condition
            sol(i,:) = yTmp;    % sufficient good solution for current step
            break;              % exit the loop
        else
            yNew = yTmp;        % new estimate for current step solution
        end
        
    end
    
    str = ['Step nr. ' num2str(i-1) ' of ' num2str(numit-1) ' done!'];
    disp(str);
    
end

%% get stationary solution
% 
%           [ lambda_11     0     ] [ r11 ]   = [ F11(t,x,R) ]
%           [    0      lambda_12 ] [ r12 ]_x = [ F12(t,x,R) ]
%

disp('Compute Stationary Solution!');
    function RHS = getRHS(R)
        if isrow(R)
            R = R';
        end
        V = ones(n, 1);
        D = diag(V, 1) - diag(V, -1);
        D(1,1:2) = [-2 2];
        D(n+1,n:n+1) = [-2 2];
        V = zeros(n+1);
        FDM = [D V V V; V D V V; V V D V; V V V D];
        
        Lambda = eye(4*n+4);
        Lambda(1:n+1,1:n+1) = diag(a*((R(1:n+1) - R(n+2:2*n+2))/2 + 1));
        Lambda(n+2:2*n+2,n+2:2*n+2) = diag(a*((R(1:n+1) - R(n+2:2*n+2))/2 - 1));
        Lambda(2*n+3:3*n+3,2*n+3:3*n+3) = diag(a*((R(2*n+3:3*n+3) - R(3*n+4:4*n+4))/2 + 1));
        Lambda(3*n+4:4*n+4,3*n+4:4*n+4) = diag(a*((R(2*n+3:3*n+3) - R(3*n+4:4*n+4))/2 - 1));
        
        FDM = Lambda * 2/deltax * FDM;
        % lower boundary
        FDM(1,1:2) = [1 0];
        FDM(1,n+2) = 1;
        % upper Boundary
        FDM(end, 3*n+3) = 1;
        FDM(end, 4*n+3:4*n+4) = [0 1];
        % coupling condition 1
        FDM(2*n+2,n+1) = -1;
        FDM(2*n+2,2*n+1:2*n+2) = [0 -1];
        FDM(2*n+2,2*n+3) = 1;
        FDM(2*n+2,3*n+4) = 1;
        % coupling condition 2
        FDM(2*n+3,n+1) = -1;
        FDM(2*n+3,2*n+2) = 1;
        FDM(2*n+3,2*n+3:2*n+4) = [control(0) 0];
        FDM(2*n+3,3*n+4) = -control(0);
        
        F = zeros(4*n+4,1);
        F(1:n+1) = -1/8*a*lambdaFric/diam * (R(1:n+1) - R(n+2:2*n+2)).^2;
        F(n+2:2*n+2) = 1/8*a*lambdaFric/diam * (R(1:n+1) - R(n+2:2*n+2)).^2;
        F(2*n+3:3*n+3) = -1/8*a*lambdaFric/diam * (R(2*n+3:3*n+3) - R(3*n+4:4*n+4)).^2;
        F(3*n+4:4*n+4) = 1/8*a*lambdaFric/diam * (R(2*n+3:3*n+3) - R(3*n+4:4*n+4)).^2;
        % lower boundary
        F(1) = 2*log(pressureLeft(0));
        % upper Boundary
        F(4*n+4) = 2*log(2*a*flowRight(0)) - 2*log(R(3*n+3) - R(4*n+4));
        % coupling condition 1
        F(2*n+2) = 2*log(control(0));
        % coupling condition 2
        F(2*n+3) = 0;
        
        RHS = F - FDM * R;
    end
    
    function DFx0 = diffquotcentralExtended(fun, x0)
        h = sqrt(eps);
        m = length(fun(x0));
        nD = length(x0); %length of the jacobimatrix
        e = eye(nD); %
        DFx0 = zeros(m, nD);
        for j = 1:nD
            DFx0(:,j) = (fun(x0 + h*e(:,j)) - fun(x0 - h*e(:,j))) / (2 * h);
        end
    end

    function [x,numit] = newtonsol(fun, x0, tol, maxit, p)
    % Initializations
    if isrow(x0)
        x = x0';
    else
        x = x0;
    end
    % get gradient
    A = diffquotcentralExtended(fun, x);
    b = fun(x);
    Anfangsresiduum = norm(fun(x0));
    grenze = tol*Anfangsresiduum;
    Residuum = Anfangsresiduum;
    numit2 = 0;
    disp('Newton Method started');
    while (grenze <= Residuum)
        numit2 = numit2 + 1;
        v = A\b; %löse Av = b
        x = x - v;
        %Jede p-Schritte wird DF neu berechnet
        if (mod(numit2, p) == 0)
            A = diffquotcentralExtended(fun, x);
        end
        b = fun(x);  
        Residuum = norm(fun(x)); 
        if numit2 == maxit
            fprintf('maxit erreicht!');
            break;
        end
        str = ['Iteration ' num2str(numit2) ' done'];
        disp(str);
    end
    str = ['Newton Method finished after ' num2str(numit2) ' iterations'];
    disp(str);
    end

fStat = @(x) getRHS(x);
solStat = newtonsol(fStat, y0, sqrt(eps), 1000, 1);
[p1Stat, p2Stat, q1Stat, q2Stat] = getPressureAndFlow(solStat);

%% plot results

    function [p1,p2,q1,q2] = getPressureAndFlow(R)
        r11 = R(1:n+1);
        r12 = R(n+2:2*n+2);
        r21 = R(2*n+3:3*n+3);
        r22 = R(3*n+4:end);
        p1 = exp((r11 + r12)/2);
        p2 = exp((r21 + r22)/2);
        q1 = (r11 - r12)/(2*a) .* exp((r11 + r12)/2);
        q2 = (r21 - r22)/(2*a) .* exp((r21 + r22)/2);
    end

figure;
xlabel('x');
ylabel('p(t,x),q(t,x)');
for i = 1:4:numit
    [p1,p2,q1,q2] = getPressureAndFlow(sol(i,:));
    subplot(3,1,1)
    plot(linspace(0,L,n+1), p1Stat, 'g', linspace(0,L,n+1), p1, 'b', linspace(L,2*L,n+1), p2Stat, 'g', linspace(L,2*L,n+1), p2, 'b');
    legend('stationary', 'dynamic');
    str = ['pressure at time t = ' num2str(tspan(i))];
    title(str);
    axis([0 2*L .9e04 2e05]);
    subplot(3,1,2)
    plot(linspace(0,L,n+1), q1Stat, 'g', linspace(0,L,n+1), q1, 'r', linspace(L,2*L,n+1), q2Stat, 'g', linspace(L,2*L,n+1), q2, 'r');
    legend('stationary', 'dynamic');
    str = ['flow at time t = ' num2str(tspan(i))];
    title(str);
    axis([0 2*L 10 50]);
    subplot(3,1,3)
    plot(linspace(0,L,n+1), a^2*q1./p1, 'm', linspace(L,2*L,n+1), a^2*q2./p2,'m');
    legend('dynamic');
    str = ['velocity at time t = ' num2str(tspan(i))];
    title(str);
    axis([0 2*L 20 60]);
    drawnow();  
end

tmp = 1;

end
