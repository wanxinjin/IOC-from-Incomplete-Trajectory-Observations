function sol=OCsolver_CostEnd(x0,horizon,dyn,path_cost,final_cost)
import casadi.*

% setup general parameter
n=dyn.size1_in('X');
m=dyn.size1_in('U');

% Extract the model parameter
controlLowerBound   =   -2000*ones(m,1);
controlUpperBound   =   2000*ones(m,1);
stateLowerBound     =   -1000*ones(n,1);
stateUpperBound     =   1000*ones(n,1);


% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
cnt={};
lbg = [];
ubg = [];



% "Lift" initial conditions
Xk = MX.sym('x0', n);
w = {w{:}, Xk};
lbw = [lbw; x0];
ubw = [ubw; x0];
w0 = [w0; x0];

% Formulate the NLP for the path
for k=1:horizon
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],m);
    w = {w{:}, Uk};
    lbw = [lbw(:); controlLowerBound];
    ubw = [ubw(:); controlUpperBound];
    w0 = [w0(:);  zeros(m,1)];
    
    % Integrate till the end of the interval
    res = dyn('X', Xk, 'U', Uk);
    Xkp = res.f;
    Ck  = path_cost('X', Xk, 'U', Uk);
    J=J+Ck.c;
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k)], n);
    w = {w{:}, Xk};
    lbw = [lbw(:);  stateLowerBound];
    ubw = [ubw(:);  stateUpperBound];
    w0 = [w0(:); zeros(n,1)];
    
    % Add equality constraint
    cnt = [cnt{:}; Xkp-Xk];
    lbg = [lbg(:);zeros(n,1)];
    ubg = [ubg(:);zeros(n,1)];
end


% adding the final cost
Ck  = final_cost('X', Xk);
J=J+Ck.c;


vertcat(w{:})
% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(cnt{:})    );
solver = nlpsol('solver', 'ipopt', prob);


% Solve the NLP
solution = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);

solution=full(solution.x);
solution=reshape(solution(n+1:end),[m+n,horizon]);

sol.u=solution(1:m,:);
sol.x=[x0,solution(m+1:end,:)];

        
end



