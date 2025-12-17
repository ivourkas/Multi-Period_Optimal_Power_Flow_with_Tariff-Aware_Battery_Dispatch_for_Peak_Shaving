function [ldl,fval,pk,btr,dspt] = BattPkShvRsCm (tm,Y,P,Q,C_sb,P_PV,ilF,P_ilM,P_invM,Q_invM,P_RChM,P_DChM,E_M,E_mt,n_b,C_if,C_PV,C_b,C_IL,ssn,yr,irr,peak)


%% ~ Input data ~
% tm:                       24 (daily planning)
% Y (n-by-n):               admittance matrix (n=number of system buses)
% P,Q (n-by-1):             maximum bus load (<0)
% C_sb (scalar OR 1-by-tm): slack bus energy cost
% P_PV (n-by-1):            installed PV capacity (25yr warrantied to down
%                           to 80% capacity)
% ilF (scalar):             factor of total load as critical - assumed
%                           flat for all
% P_ilM (n-by-1):           maximum interrupted bus load for scheduling 
%                           horizon ahead
% P_invM (n-by-1):          inverter rated/maximum active power
% Q_invM (n-by-1):          inverter rated/maximum reactive power
% P_RChM (n-by-1):          charger maximum absorbed active power 
% P_DChM (n-by-1):          charger maximum injected active power
% E_M (n-by-1):             battery capacity (precalculated w/ degradation)
% E_mt (n-by-1):            level of battery charging prior to scheduling 
%                           horizon
% n_b (n-by-1):             battery efficiency
% C_if (n-by-tm OR scalar): infeed cost of energy PER BUS PER HOUR (MODIFIED)
% C_PV (n-by-1 OR scalar): contracted PV energy price per bus (MODIFIED)
% C_b (scalar):             battery degradation cost (precalculated w/ 
%                           degradation)
% C_IL (n-by-tm OR 1-by-tm): interruptible load value per bus per hour
% ssn (scalar):             season (1: summer/spring, 2: winter/fall)
% yr  (scalar >=1):         year from the start of investment
% irr (yr-by-1):            investment return rate or inflation
% peak:                     pareto front gain to assess quasi peak 
%                           shaving effect


%% ~ Vector of Variables in the Linear Program ~
% n buses:                  "1" slack bus, "N" last bus
% x1:                       vector of variables for hour t=1
% x1=[P1,P2i,P2ii,P2iii,P2iv,P2v,P2vi,P3i,...,PNv,PNvi,d1,d2,...,dN,P1.2,
% P1.3,P2.3,...,PN-1.N,P2.1,..,PN.N-1,..,Q1,Q2i,Q2ii,Q2iii,Q2iv,Q2v,Q2vi,
% Q3i,...,QNv,QNvi,V1,V2,...,Vn,Q1.2,Q1.3,Q2.3,...,QN-1.N,Q2.1,..,QN.N-1]'
% concatenated for t=[1,horizon ahead end] 
% i.e. [slack bus active power, [bus active power from actors (i,ii,...), 
% bus load] concatenated for n buses, [volt angle] for n buses, [active 
% flows] for n buses, slack bus reactive power, [bus reactive power from  
% actors (i,ii,...), bus load] for n buses, [volt magnitude] for n buses,  
% [active flows] for n buses,]'
% actors break-down: (i) net-metering PV, (ii) FiT PV, (iii) Battery 
% discharge, (iv) Battery Recharge [<0], (v) interruptible load


%% ~ Equality Constraints Matrix ~
% Aeq1:                     equality constraints matrix for hour t=1
% Aeq1=[active power flows to buses; active power flows constraints; 
% reactive power flows to buses; mixed power flows constraints;]
% concatenated for t=[1,horizon ahead end] and then append (n-1)*t lines 
% equality constraints for reactive power of interrupted load 

% 
%% ~ Inequality Constraints Matrix ~
% Aineq=[[PV sum];[horizon IL];[battery capacity];[inverter limits];[quasi
% peak shaving]]
% PV sum ((n-1)*t lines):   scheduled PV on net-metering and contracted
%                           price not more than available
% horizon IL (n-1 lines):   total interruptible load in the horizon
% battery capacity (2*(n-1)*t lines):   re-/dis-charging added/subtracted
%                                       to previous horizon up to total
% inverter limits (2*(n-1)*t lines):    min-max of active & reactive power
% quasi peak shaving (2*t-1 lines):     valley filling by minimizing 
%                                       deviation from load average

%% initializations
INF=1e6;
n=length(Y);

%% INPUT VALIDATION - CHECK BATTERY STATES
for j = 2:n
    if E_mt(j) > E_M(j)
        fprintf('WARNING: Bus %d has E_mt(%.6f) > E_M(%.6f) at function entry!\n', j, E_mt(j), E_M(j));
        E_mt(j) = E_M(j);
    end
    if E_M(j) < 0 || E_mt(j) < 0
        fprintf('WARNING: Bus %d has negative battery values: E_M=%.6f, E_mt=%.6f\n', j, E_M(j), E_mt(j));
        E_M(j) = max(0, E_M(j));
        E_mt(j) = max(0, E_mt(j));
    end
end

P=-P; Q=-Q; 
% change horizon 
Pl=zeros(n,1,tm); 
cl=Pl;
P_lM=P; acosf=sign(P).*sign(Q).*acos(abs(P)./(sqrt(P.^2+Q.^2))); 
% load initialize

PgU=zeros(n,5,tm); PgD=PgU; QgU=PgU; QgD=PgU; cg=PgU;
% cg:   n rows from nodes // 5 columns for prices of (i) net-metering PV, 
%       (ii) FiT PV, (iii) Battery discharge, (iv) Battery Recharge [<0], 
%       (v) interruptible load // tm aisles for planning horizon

% Pb=zeros(n*tm*5,1);
% Pbd=Pb;

% Compound inflation to current year ONCE OUTSIDE THE LOOP!
C_b_inflated = C_b * prod(1 + irr(1:yr)); 
% Handle C_if - make it per-bus if not already
if isscalar(C_if)
    C_if = repmat(C_if, n, tm);  % Expand scalar to matrix
elseif size(C_if,1) == 1 && size(C_if,2) == tm
    % 1×24 vector - expand to all buses
    C_if = repmat(C_if, n, 1);
end
% Now C_if is n×tm matrix

% Handle C_PV - make it per-bus if not already
if isscalar(C_PV)
    C_PV = C_PV * ones(n, 1);  % Expand scalar to vector
end
% Now C_PV is n×1 vector

if isscalar(C_IL)
    C_IL = repmat(C_IL, n, tm);  % Expand scalar to n×tm matrix
elseif size(C_IL,1) == 1 && size(C_IL,2) == tm
    % 1×24 vector - expand to all buses (backward compatibility)
    C_IL = repmat(C_IL, n, 1);
end
% Now C_IL is n×tm



% change start time
startt=round(rand*24, 0);
for (i=1:1:tm)
    exi=startt+i;
    if exi>24
        exi=-(24-startt-i);
    end
    for (j=1:1:n)
        buf=(0.8^(yr/25))*(1-0.5*round(rand*0.3*ssn))*min(1,max(0,min(exi-1-6,18-exi+1)))*max(0,(max(1,max(0,min(exi-1-6,18-exi+1)))/6-(0.2+(ssn-1)*0.6)*rand));
        PgU(j,1,i)=buf*P_PV(j);
        PgU(j,2,i)=buf*P_PV(j);
        if (buf>0)
            QgU(j,1,i)=Q_invM(j);
            QgU(j,2,i)=Q_invM(j);
            QgD(j,1,i)=-Q_invM(j);
            QgD(j,2,i)=-Q_invM(j);
        end
        % PV (NM & Cntrl)
        hr=exi-1;
        buf=(-1.6*sin(hr*pi/10)+2.4)*min(1,max(0,min(hr+1,11-hr)))*(0.5+0.2*rand)+(-0.8*sin(-8*pi/6+hr*pi/7)+3.2)*min(1,max(0,min(hr-10,24-hr)))*(0.6+0.25*rand);
        Pl(j,1,i)=buf*P_lM(j)/4;
        %if (i>1)
        %    if (lb(copf_vars+(i-1)*10*n+3*n+j,1)<lb(copf_vars+(i-2)*10*n+3*n+j,1))
        %        for (k=0:1:i-2)
        %            lb(copf_vars+k*10*n+3*n+j,1)=lb(copf_vars+(i-1)*10*n+3*n+j,1);
        %        end
        %    end
        %end 
        % minimum critical of all hours across all hours ? ? ?
        % ~~~ 
        % load
        PgU(j,3,i)=min(P_invM(j),P_DChM(j));
        PgD(j,4,i)=-min(P_invM(j),P_RChM(j));
        QgU(j,3,i)=Q_invM(j);
        QgU(j,4,i)=Q_invM(j);
        QgD(j,3,i)=-Q_invM(j);
        QgD(j,4,i)=-Q_invM(j);
        % cg(j,4,i) = C_if(j,i) + C_b_inflated;  % Use C_if(j,i) instead of C_if(i)
        
        % DCh & RCh
        PgU(j,5,i)=Pl(j,1,i)-min(ilF(j)*P_lM(j)/4,Pl(j,1,i));
        % instead of PgU(j,5,i)=Pl(j,1,i)-min(ilF*P_lM(j)/4,Pl(j,1,i));
        QgU(j,5,i)=Pl(j,1,i);
        QgD(j,5,i)=-Pl(j,1,i);
        % % IL 
        % % Feed-in tariff - NOW PER-BUS
        % for (j=1:1:n)
        %     cg(j,2,:) = -C_PV(j);  % Use C_PV(j) instead of scalar C_PV
        % end

        % % Net-metering - uses C_if per-bus
        % for (i=1:1:tm)
        %     for (j=1:1:n)
        %         cg(j,1,i) = -C_if(j,i);  % Use C_if(j,i)
        %         cg(j,3,i) = -C_if(j,i);  % Battery discharge value
        %     end
        % end
        % cg(j,5,:)=-C_IL;    % Negative = revenue from shedding  (%C_IL-C_if;)
        cl(j,:,:)=0; %C_if; no-need; already cost on Slack Bus
        % add comments
    end
end
%% SET COSTS - AFTER the main loop, with correct indexing
% Battery recharge cost (per-bus, per-hour)
for i=1:tm
    for j=1:n
        cg(j,4,i) = -(C_if(j,i) + C_b_inflated);
    end
end

% Feed-in tariff (per-bus, all hours)
for j=1:n
    for i=1:tm
        cg(j,2,i) = -C_PV(j);
    end
end

% Net-metering and battery discharge (per-bus, per-hour)
for i=1:tm
    for j=1:n
        cg(j,1,i) = -C_if(j,i);  % Net-metering
        cg(j,3,i) = -C_if(j,i);  % Battery discharge value
    end
end

% Interruptible load
for i=1:tm
    for j=1:n
        cg(j,5,i) = -C_IL(j,i);  
    end
end

% Slack bus cost (overwrites generator 1 cost for bus 1)
for i=1:tm
    cg(1,1,i) = C_sb(i);
end

%% =========================================================================
%% SOFT CONSTRAINT: Discourage simultaneous charge/discharge
%% =========================================================================
% Instead of MILP binary variables, we add a small cost penalty that makes
% simultaneous operation economically suboptimal.
%
% The penalty must be:
% - Large enough to prevent simultaneous operation
% - Small enough to not distort the primary optimization objective
%
% We use a penalty proportional to the efficiency loss that would occur

% Calculate penalty based on round-trip efficiency loss
% If battery charges and discharges simultaneously, energy is wasted
% Penalty = cost of wasted energy

simultaneous_penalty_factor = 0.01;  % 1% of grid price as penalty

for i = 1:tm
    for j = 2:n
        if E_M(j) > 0  % Only for buses with batteries
            % Add small penalty to discharge cost
            % This breaks ties when optimizer might consider simultaneous operation
            current_grid_price = mean(C_if(j, :));  % Average grid price for this bus
            penalty = simultaneous_penalty_factor * current_grid_price;
            
            % Discharge cost adjustment (make it slightly more expensive)
            % cg(j,3,i) is discharge credit (negative), so we make it less negative
            cg(j, 3, i) = cg(j, 3, i) + penalty;
            
            % Charge cost adjustment (make it slightly more expensive)  
            % cg(j,4,i) is charge cost (negative for expense, recharge output negative), so we make it more negative
            cg(j, 4, i) = cg(j, 4, i) - penalty;
        end
    end
end

QgU(2:end,:,:)=0; % zeroed all except slack bus (keep?)
QgD(2:end,:,:)=0; % zeroed all except slack bus (keep?)
% costs, ubs & lbs...


%% OPF modeling
[Aeq,beq,lb,ub,ct,lp_vars] = COPF_fe(Y,Pl,acosf,cl,PgU,PgD,QgU,QgD,cg);

% Display ct information
% disp('Size of ct:');
% disp(size(ct));

% Explanation given from COPF_fe:
% Buses: n = 10 (from Y matrix)
% Components per Bus (non-slack):
% ->Active Power: 5 (PV net-metering, PV FiT, Battery discharge, Battery recharge, Interruptible load).
% ->Fixed Load: 1 (non-dispatchable).
% Lines: s_flows = 9 (s_flows=length(nonzeros(triu(Y,1)));).
% Thus, lp_vars = 1 + (n-1)*(wPg + wPl) + n + 2*s_flows;
% -> 1: Slack bus active power (P1).
% -> (n-1)*(wPg + wPl): Non-slack buses’ active components (wPg = 5 for PV/battery/IL, wPl = 1 for load).
% -> n: Voltage angles (d1..d10).
% -> 2*s_flows: Active power flows (forward and reverse for each line).
%
% Finally, for n = 10, s_flows = 9:
% 
% lp_vars = 1 + 9*(5+1) + 10 + 2*9 = 1 + 54 + 10 + 18 = 83.
% Total Variables for 24 Hours
% -> Active Power Variables/Hour: 83.
% -> Reactive Power Variables/Hour: 83 (structured similarly to active power).
% 
% Total Variables/Hour: 83 (active) + 83 (reactive) = 166.
% 
% Total Variables for 24h: 166 × 24 = 3984. CONFIMED


%% IL reactive according to cosf(j)
if (sum(sum(sum(abs(QgU(2:end,:,:)))))>0)
    nxAeq=size(Aeq,1)+1
    for (i=1:1:tm)
        for (j=2:1:n)
            Aeq(nxAeq,(i-1)*lp_vars*2+lp_vars+1+(j-2)*6+5)=1;
            Aeq(nxAeq,(i-1)*lp_vars*2+1+(j-2)*6+5)=-tan(acosf(j));
            beq(nxAeq,1)=0;
            nxAeq=nxAeq+1;
        end
    end
end


%% PV Net-metering & PV FiT not greater than available
Aineq=zeros(1,tm*lp_vars*2);
bineq = []; 
for (i=1:1:tm)
    for (j=2:1:n)
        Aineq((i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+1)=1;
        Aineq((i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+2)=1;
        bineq((i-1)*(n-1)+j-1,1)=ub((i-1)*lp_vars*2+1+(j-2)*6+1,1);
    end
end
for i=size(bineq,1):-1:1
    if ~(abs(bineq(i,:))>0)
        Aineq(i,:)=[];
        bineq(i,:)=[];
    end
end


%% max interrupted/responsive load demand per day
Aineq1=zeros(1,tm*lp_vars*2);
bineq1=[];
for (j=2:1:n)
    Aineq1(j-1,1+(j-2)*6+5:lp_vars*2:end)=1;
    bineq1(j-1,1)=P_ilM(j);
    % max shed IL per day
end
for i=size(bineq1,1):-1:1
    if ~(abs(bineq1(i,:))>0)
        Aineq1(i,:)=[];
        bineq1(i,:)=[];
    end
end
Aineq=[Aineq;Aineq1];
bineq=[bineq;bineq1];


%% battery discharging-charging less than E_M-E_mt (right side constraint)
% Energy balance: E(t) = E(t-1) - discharge/η + charge*η
% Constraint: cumulative (charge*η - discharge/η) ≤ E_M - E_mt
PV_load_ineq1=size(Aineq,1);
Aineq1=zeros(1,tm*lp_vars*2);
bineq1=[];
for (i=1:1:tm)
    for (j=2:1:n)
        for (k=1:1:i)
            % Discharge removes energy from battery (divided by efficiency)
            Aineq1((i-1)*(n-1)+j-1,(k-1)*lp_vars*2+1+(j-2)*6+3) = -1/n_b(j);
            % Charge adds energy to battery (multiplied by efficiency)
            Aineq1((i-1)*(n-1)+j-1,(k-1)*lp_vars*2+1+(j-2)*6+4) = -1*n_b(j);
        end
        bineq1((i-1)*(n-1)+j-1,1) = E_M(j) - E_mt(j);
    end
end
bineq2=repmat(E_mt(2:n,1),tm,1);
for i=size(bineq1,1):-1:1
    if ~(abs(bineq1(i,:))>0)
        Aineq1(i,:)=[];
        bineq1(i,:)=[];
        bineq2(i,:)=[];
    end
end
Aineq=[Aineq;Aineq1];
bineq=[bineq;bineq1];

%% battery discharging-charging greater than E_mt (left side constraint)
PV_load_batt_ineq=size(Aineq,1); 
Aineq=[Aineq;-Aineq(PV_load_ineq1+1:PV_load_batt_ineq,:)];
bineq=[bineq;bineq2];

%% =========================================================================
%% MINIMUM SOC CONSTRAINT - Battery cannot discharge below SOC_min
%% =========================================================================
% Physical meaning: E(t) >= SOC_min * E_M for all time periods
% This prevents deep discharge which damages battery life
%
% Energy at time t: E(t) = E_mt + sum_{k=1}^{t}(charge_k * n_b - discharge_k / n_b)
% Constraint: E(t) >= SOC_min * E_M
% Rearranged: sum_{k=1}^{t}(discharge_k / n_b - charge_k * n_b) <= E_mt - SOC_min * E_M

SOC_min = 0.10;  % 10% minimum state of charge - ADJUSTABLE PARAMETER

% Count existing inequality constraints
PV_load_batt_ineq_SOC_start = size(Aineq, 1);

% Build SOC floor constraints for each bus and each hour
Aineq_SOC = zeros((n-1)*tm, size(Aineq, 2));
bineq_SOC = zeros((n-1)*tm, 1);

constraint_row = 0;
for i = 1:tm  % For each hour
    for j = 2:n  % For each non-slack bus
        constraint_row = constraint_row + 1;
        
        if E_M(j) > 0  % Only for buses with batteries
            % Cumulative constraint: sum of (discharge/η - charge*η) <= E_mt - SOC_min*E_M
            for k = 1:i  % Sum from hour 1 to current hour i
                % Discharge variable index (generator 3)
                discharge_idx = (k-1)*lp_vars*2 + 1 + (j-2)*6 + 3;
                % Charge variable index (generator 4)
                charge_idx = (k-1)*lp_vars*2 + 1 + (j-2)*6 + 4;
                
                % Discharge drains battery: positive coefficient (uses energy)
                % Note: discharge is positive in the model
                Aineq_SOC(constraint_row, discharge_idx) = 1 / n_b(j);
                
                % Charge fills battery: negative coefficient (adds energy)
                % Note: charge is negative in the model, so we use positive coeff
                % to make it: -(-charge)*η = +charge*η subtracted from drain
                Aineq_SOC(constraint_row, charge_idx) = 1 * n_b(j);
            end
            
            % Right-hand side: maximum allowed drain = E_mt - SOC_min * E_M
            bineq_SOC(constraint_row) = E_mt(j) - SOC_min * E_M(j);
        else
            % No battery: constraint is trivially satisfied (0 <= 0)
            bineq_SOC(constraint_row) = 0;
        end
    end
end

% Remove rows for buses with no battery (all zeros)
rows_to_remove = [];
for i = size(bineq_SOC, 1):-1:1
    if all(Aineq_SOC(i, :) == 0) && bineq_SOC(i) == 0
        rows_to_remove = [rows_to_remove, i];
    end
end
Aineq_SOC(rows_to_remove, :) = [];
bineq_SOC(rows_to_remove) = [];

% Append to main inequality constraints
if ~isempty(bineq_SOC)
    Aineq = [Aineq; Aineq_SOC];
    bineq = [bineq; bineq_SOC];
end

% Display constraint info (can be commented out for production)
% fprintf('Added %d minimum SOC constraints (SOC_min = %.0f%%)\n', size(Aineq_SOC, 1), SOC_min*100);

%% net-metering, FiT, battery discharging+recharging not greater than 
% inverter max active power & within +/- max of reactive
PV_load_batt_ineq=size(Aineq,1);
for (i=1:1:tm)
    for (j=2:1:n)
        Aineq(PV_load_batt_ineq+(i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+1)=1;
        Aineq(PV_load_batt_ineq+(i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+2)=1;
        Aineq(PV_load_batt_ineq+(i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+3)=1;
        Aineq(PV_load_batt_ineq+(i-1)*(n-1)+j-1,(i-1)*lp_vars*2+1+(j-2)*6+4)=1;
        bineq(PV_load_batt_ineq+(i-1)*(n-1)+j-1,1)=P_invM(j);
        Aineq(PV_load_batt_ineq+(i-1+tm)*(n-1)+j-1,(i-1)*lp_vars*2+lp_vars+1+(j-2)*6+1)=1;
        Aineq(PV_load_batt_ineq+(i-1+tm)*(n-1)+j-1,(i-1)*lp_vars*2+lp_vars+1+(j-2)*6+2)=1;
        Aineq(PV_load_batt_ineq+(i-1+tm)*(n-1)+j-1,(i-1)*lp_vars*2+lp_vars+1+(j-2)*6+3)=1;
        Aineq(PV_load_batt_ineq+(i-1+tm)*(n-1)+j-1,(i-1)*lp_vars*2+lp_vars+1+(j-2)*6+4)=1;
        bineq(PV_load_batt_ineq+(i-1+tm)*(n-1)+j-1,1)=Q_invM(j);
    end
end
Aineq=[Aineq;-Aineq(PV_load_batt_ineq+tm*(n-1)+1:PV_load_batt_ineq+2*tm*(n-1),:)];
bineq=[bineq;bineq(PV_load_batt_ineq+tm*(n-1)+1:PV_load_batt_ineq+2*tm*(n-1),1)];
for i=size(bineq,1):-1:PV_load_batt_ineq+1
    if ~(abs(bineq(i,:))>0)
        Aineq(i,:)=[];
        bineq(i,:)=[];
    end
end
%size(Aineq)
%size(Aeq)


%% quasi-peak-shaving
PV_load_batt_inv_ineq=size(Aineq,1);
aux_ineq=size(Aineq,2);
for (i=1:1:tm-1)
    Aineq(PV_load_batt_inv_ineq+2*i-1,(i-1)*lp_vars*2+1)=1;
    Aineq(PV_load_batt_inv_ineq+2*i-1,i*lp_vars*2+1)=-1;
    Aineq(PV_load_batt_inv_ineq+2*i-1,aux_ineq+i)=-1;
    bineq(PV_load_batt_inv_ineq+2*i-1,1)=0;
    Aineq(PV_load_batt_inv_ineq+2*i,(i-1)*lp_vars*2+1)=-1;
    Aineq(PV_load_batt_inv_ineq+2*i,i*lp_vars*2+1)=1;
    Aineq(PV_load_batt_inv_ineq+2*i,aux_ineq+i)=-1;
    bineq(PV_load_batt_inv_ineq+2*i,1)=0;
end
Aeq=[Aeq,zeros(size(Aeq,1),tm-1)];
lb=[lb;zeros(tm-1,1)];
ub=[ub;1000*ones(tm-1,1)];
ct=[ct;peak*ones(tm-1,1)];
%size(Aineq)
%size(Aeq)
% 
% format short 
% idx=(1:length(lb))';
% fprintf('Index \t Lower Bound \t Upper Bound\n');
% fprintf('----- \t ----------- \t -----------\n');
% for i = 1:length(lb)
%     fprintf('%5d \t %10.4f \t %10.4f\n', idx(i), lb(i), ub(i));
% end

%% CRITICAL FIX: Force battery variable bounds to zero when capacity is zero
for j = 2:n  % Skip slack bus
    if E_M(j) == 0  % If battery capacity is zero
        for i = 1:tm  % For each time step
            % Calculate variable indices for battery discharge and charge
            discharge_var_idx = (i-1)*lp_vars*2 + 1 + (j-2)*6 + 3;  % Generator 3
            charge_var_idx = (i-1)*lp_vars*2 + 1 + (j-2)*6 + 4;     % Generator 4
            
            % Force both lower and upper bounds to exactly zero
            lb(discharge_var_idx) = 0;
            ub(discharge_var_idx) = 0;
            lb(charge_var_idx) = 0;
            ub(charge_var_idx) = 0;
        end
    end
end

% options = optimoptions('linprog','Algorithm','dual-simplex');%,'Display','none');
% [x,fval]=linprog(ct,Aineq,bineq,Aeq,beq,lb,ub,options);%,[]); 
%% CRITICAL FIX: Add robust error handling for failed/incomplete solutions
options = optimoptions('linprog', ...
    'Algorithm', 'dual-simplex', ...
    'Display', 'off', ...
    'ConstraintTolerance', 1e-6, ...      % Looser tolerance
    'OptimalityTolerance', 1e-6, ...
    'MaxIterations', 10000); 
% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[x, fval, exitflag, output] = linprog(ct, Aineq, bineq, Aeq, beq, lb, ub, options);

% ERROR HANDLING - CHECK IF SOLVER FAILED OR RETURNED INCOMPLETE SOLUTION
if isempty(x) || exitflag <= 0 || length(x) ~= (2*lp_vars*tm + (tm-1))
    fprintf('\n========== LINPROG FAILED OR INCOMPLETE ==========\n');
    fprintf('Exit flag: %d\n', exitflag);
    if ~isempty(output)
        fprintf('Message: %s\n', output.message);
    end
    fprintf('Expected solution size: %d\n', 2*lp_vars*tm + (tm-1));
    fprintf('Actual solution size: %d\n', length(x));
    fprintf('Year: %d, Season: %d\n', yr, ssn);
    
    % Check problematic battery states
    fprintf('Battery states at failure:\n');
    for j = 2:n
        if E_M(j) > 0
            fprintf('  Bus %d: E_M=%.6f, E_mt=%.6f\n', j, E_M(j), E_mt(j));
        end
    end
    fprintf('====================================\n\n');
    
    % Return safe default values instead of crashing
    ldl = 0;
    fval = Inf;  % High cost to indicate infeasibility
    pk = zeros(n,1);
    btr = 0;
    dspt = zeros(1 + (n-1)*5, tm);  % Default dispatch matrix
    
    % Early return to avoid processing invalid solution
    return;
end
% options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
% [x,fval,exitflag,output]=linprog(ct,Aineq,bineq,Aeq,beq,lb,ub,options);
% 

% % Check solution vector size
% expected_size = 2*lp_vars*tm + (tm-1);  % tm-1 for peak shaving auxiliary variables
% if length(x) ~= expected_size
%     fprintf('\n========== SOLUTION SIZE MISMATCH ==========\n');
%     fprintf('Expected size: %d, Got: %d\n', expected_size, length(x));
%     fprintf('Year: %d, Season: %d\n', yr, ssn);
%     fprintf('==========================================\n\n');
% 
%     ldl = NaN;
%     fval = Inf;
%     pk = zeros(n-1,1);
%     btr = NaN;
%     dspt = [];
%     return;
% end

% DAY-AHEAD SCHEDULING OUTPUT
% [x,fval]=linprog(ct,[],[],Aeq,beq,lb,ub,[],options);
% DAY-AHEAD SCHEDULING OUTPUT





%% --- Calculate Power Balance ---

% 1. Print the solution vector x with indices in two columns
% Export just the values (single column)
% filename_values = 'x_values.xlsx';
% writematrix(x, filename_values);
% fprintf('Values exported to %s\n', filename_values);

% % 2. Verify power balance for each hour
% n_hours = 24;
% vars_per_hour = 166;
% 
% for h = 1:n_hours
%     hour_offset = (h-1) * vars_per_hour;
% 
%     % Active power balance check (components 1-55 in the hour's block)
%     active_start = hour_offset + 1;
%     active_end = hour_offset + 55;
%     active_sum = sum(x(active_start:active_end));
% 
%     % Reactive power balance check (components 184-238 in full vector = 83+1:83+55 in hour's block)
%     reactive_start = hour_offset + 83 + 1;
%     reactive_end = hour_offset + 83 + 55;
%     reactive_sum = sum(x(reactive_start:reactive_end));
% 
%     fprintf('\nHour %d Power Balance:\n', h);
%     fprintf('Active power sum (components %d-%d): %f\n', active_start, active_end, active_sum);
%     fprintf('Reactive power sum (components %d-%d): %f\n', reactive_start, reactive_end, reactive_sum);
% 
%     % Check if balances are zero (within numerical tolerance)
%     tol = 1e-6;
%     if abs(active_sum) > tol
%         fprintf('Warning: Active power imbalance detected in hour %d\n', h);
%     end
%     if abs(reactive_sum) > tol
%         fprintf('Warning: Reactive power imbalance detected in hour %d\n', h);
%     end
% end



% %% 3. Battery Power Balance Verification
% n_hours = 24;
% vars_per_hour = 166;
% non_slack_buses = 2:10; % Exclude slack bus (bus 1)
% tol = 1e-6; % Numerical tolerance
% 
% % Component positions (0-based index per bus)
% BATT_CHARGE_P = 3;   % 3rd component = Battery charging (P)
% BATT_DISCHARGE_P = 4; % 4th component = Battery discharging (P)
% BATT_CHARGE_Q = BATT_CHARGE_P + 83;  % Reactive power charging
% BATT_DISCHARGE_Q = BATT_DISCHARGE_P + 83; % Reactive power discharging
% 
% fprintf('\nBattery Power Balance Verification (Charge + Discharge = 0):\n');
% fprintf('%-6s %-12s %-12s %-12s %-12s\n', ...
%     'Hour', 'P_Sum', 'Balanced?', 'Q_Sum', 'Balanced?');
% 
% for h = 1:n_hours
%     hour_offset = (h-1) * vars_per_hour;
%     P_sum = 0;
%     Q_sum = 0;
% 
%     % Calculate net battery power across all non-slack buses
%     for k = non_slack_buses
%         % Active power sum (Charge + Discharge)
%         P_sum = P_sum + x(hour_offset + 1 + (k-2)*6 + BATT_CHARGE_P) + ...
%                        x(hour_offset + 1 + (k-2)*6 + BATT_DISCHARGE_P);
% 
%         % Reactive power sum (Charge + Discharge)
%         Q_sum = Q_sum + x(hour_offset + 1 + (k-2)*6 + BATT_CHARGE_Q) + ...
%                        x(hour_offset + 1 + (k-2)*6 + BATT_DISCHARGE_Q);
%     end
% 
%     % Check balance
%     P_balanced = abs(P_sum) < tol;
%     Q_balanced = abs(Q_sum) < tol;
% 
%     % Print results
%     fprintf('%-6d %-12.6f %-12s %-12.6f %-12s\n', ...
%         h, P_sum, string(P_balanced), Q_sum, string(Q_balanced));
% 
%     % Detailed warnings if imbalance exists
%     if ~P_balanced
%         fprintf('  -> Active power imbalance: %.2e\n', P_sum);
%     end
%     if ~Q_balanced
%         fprintf('  -> Reactive power imbalance: %.2e\n', Q_sum);
%     end
% end

%% x structure: tm times concatenated of
% [Slack bus 1 P, [bus 2 P of PV net-metering, PV FiT, Batt DCh, Batt RCh,
% DR, load], [bus 3 P...], ..., [bus n P...], Slack bus 1 Q, [bus 2 Q...],
% ..., [bus n Q...]]


il=[]; ldrP=[]; btP=[]; genP=[];
for (k=1:1:tm)
    for (j=1:1:5)
        for (i=1:1:n)
            if (i>1)
                if (j==5)
                    il(1+(i-2)*5+j,k)=x((k-1)*lp_vars*2+1+(i-2)*6+j,1);
                elseif (j==3)
                    btP(1+(i-2)*5+j,k)=x((k-1)*lp_vars*2+1+(i-2)*6+j,1);
                end
                ldrP(i-1,k)=x((k-1)*lp_vars*2+1+(i-1)*6,1);
                genP(1+(i-2)*5+j,k)=x((k-1)*lp_vars*2+1+(i-2)*6+j,1);
                ldrP(i-1,k)=x((k-1)*lp_vars*2+1+(i-1)*6,1);
                genQ(1+(i-2)*5+j,k)=x((k-1)*lp_vars*2+lp_vars+1+(i-2)*6+j,1);
                ldrQ(i-1,k)=x((k-1)*lp_vars*2+lp_vars+1+(i-1)*6,1);
            else
                genP(1,k)=x((k-1)*lp_vars*2+1,1);
                genQ(1,k)=x((k-1)*lp_vars*2+lp_vars+1,1);
            end
        end
    end
     voltr(:,2*(k-1)+1:2*k)=[x((k-1)*2*lp_vars+lp_vars+1+(n-1)*6+1:(k-1)*2*lp_vars+lp_vars+1+(n-1)*6+n),x((k-1)*2*lp_vars+1+(n-1)*6+1:(k-1)*2*lp_vars+1+(n-1)*6+n)];
end
ldl=sum(ldrP)-sum(il);
ldl=sum(ldl);
btr=sum(btP);
btr=sum(btr);
% fval
% x(size(x,1)-12+1:size(x,1),1)
fval=fval-peak*sum(x(size(x,1)-tm+1+1:size(x,1),1));
% grid energy cost by subtracting peak-shaving pseudo-cost 
% ! REMOVED and cost of slack bus energy (-sum(C_sb.*genP(1,:))) REMOVED !

dspt=genP;

% Initialize a matrix to store the peak calculation for each bus and each hour
% Rows represent buses (1 to n), columns represent time steps (1 to tm)
pk_bus_hourly_values = zeros(n, tm);

% Loop through each hour (time step) of the planning horizon
for k = 1:tm
    % SLACK BUS (Bus 1) - Peak is just its own load
    % For the slack bus, the peak is simply the power it supplies to meet system demand
    pk_bus_hourly_values(1, k) = abs(genP(1, k)); % Slack bus active power
    
    % NON-SLACK BUSES (Bus 2 to n) - Calculate individual bus peaks
    for i_bus = 2:n % 'i_bus' represents the actual bus index
        
        % Get the bus load for this specific bus
        bus_load_at_hour_k = ldrP(i_bus-1, k); % ldrP is indexed from 1 for bus 2
        
        % Determine the base row index in the 'genP' matrix for the current bus's generators
        base_genP_row_for_current_bus = 1 + (i_bus-2)*5;
        
        % Get individual generator contributions for this specific bus
        bus_recharge_batt = abs(genP(base_genP_row_for_current_bus + 4, k)); % j=4: Battery recharge (negative, so take abs)
        bus_net_met_PV = genP(base_genP_row_for_current_bus + 1, k);          % j=1: Net-metering PV
        bus_feed_in_tariff_PV = genP(base_genP_row_for_current_bus + 2, k);   % j=2: Feed-in-Tariff PV
        bus_discharge_batt = genP(base_genP_row_for_current_bus + 3, k);       % j=3: Battery Discharge
        bus_interr_load = genP(base_genP_row_for_current_bus + 5, k);          % j=5: Interruptible Load
        
        % Calculate the 'pk' value for this specific bus at hour 'k' using the formula:
        % pk_bus = |bus load| + |recharge batt| - |net met PV| - |feed in tariff PV| - |discharge batt| - |interr. load|
        pk_bus_hourly_values(i_bus, k) = bus_load_at_hour_k + bus_recharge_batt - ...
                                         bus_net_met_PV - bus_feed_in_tariff_PV - ...
                                         bus_discharge_batt - bus_interr_load;
    end
end

% Calculate the peak for each bus (maximum value across all hours for each bus)
pk = max(pk_bus_hourly_values, [], 2); % Max along dimension 2 (across hours)


% % Sanity check for battery energy balance
% total_battery_usage = 0;
% 
% for i = 2:10  % For each bus from 2 to 10
%     for k = 1:24  % For each hour
%         % Generator 3 is battery discharging (positive when discharging)
%         battery_net =  dspt((i-2)*5 + 4, k);
%         total_battery_usage = total_battery_usage + battery_net;
%     end
% end
% 
% % Compare with calculated btr
% disp(['Sum of battery discharging across all buses and hours: ', num2str(total_battery_usage)]);
% disp(['Calculated total battery energy used (btr): ', num2str(btr)]);
% 
% if abs(total_battery_usage - btr) < 1e-6  % Allow for small numerical differences
%     disp('Sanity check PASSED: Battery energy balances match');
% else
%     disp('Sanity check FAILED: Battery energy balances DO NOT match');
%     disp(['Difference: ', num2str(total_battery_usage - btr)]);
% end

% s=strcat('Aeqlt',datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([Aeq,beq]),'-append','delimiter','\t','newline','pc','precision','%.7f');
% s=strcat('Aineqlt',datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([Aineq,bineq]),'-append','delimiter','\t','newline','pc','precision','%.7f');
% s=strcat('ublbct',datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([ub';lb']),'-append','delimiter','\t','newline','pc','precision','%.7f');

% s=strcat('V_sb',num2str(C_sb),'PV',num2str(sum(P_PV)/sum(P)),datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat(voltr),'-append','delimiter','\t','newline','pc','precision','%.7f');
% s=strcat('P_sb',num2str(C_sb),'PV',num2str(sum(P_PV)/sum(P)),datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([genP;[];ldrP;[];fval*ones(1,tm)]),'-append','delimiter','\t','newline','pc','precision','%.7f');%
% s=strcat('Q_sb',num2str(C_sb),'PV',num2str(sum(P_PV)/sum(P)),datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([genQ;ldrQ]),'-append','delimiter','\t','newline','pc','precision','%.7f');
% s=strcat('ct_sb',num2str(C_sb),'PV',num2str(sum(P_PV)/sum(P)),datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat(ct),'-append','delimiter','\t','newline','pc','precision','%.7f');


end