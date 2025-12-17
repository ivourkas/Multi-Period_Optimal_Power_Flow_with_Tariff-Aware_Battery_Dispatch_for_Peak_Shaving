    function [Aeq,beq,lb,ub,ct,lp_vars] = COPF_fe(Y,Pl,acosf,cl,PgU,PgD,QgU,QgD,cg)

% P,Q	: active, reactive power
% ~~~~~~~~~ 3d input ~~~~~~~~~
% ~ x-d : bus
% ~ y-d : component on bus
% ~ z-d : time slot
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% U,D   : upper, lower limmits
% l,g   : load, generator
% cosf  : load power factor

% buses: 0,1,2,...,n
% x=[P1,[P2],...,[Pn],d1,...,dn,P0.1,...,Pn-1.n,..,P1.0,..,Pn.n-1,
% Q1,[Q2],...,[Qn],V1,...,Vn,Q0.1,...,Qn-1.n,Q1.0,...,Qn.n-1]'
% i.e. [bus active power, bus angle, active flows (vertical scan triu(Y)),
% bus reactive power, bus vltgs, reactive flows (vertical scan triu(Y))]'
% IMPORTANT buses power = [slack bus 1 load/gen, [bus 2 net-metering PV, 
% FiT PV, Battery discharge, Battery Recharge, interruptible load, load],
% [bus 3...], ...]

% Aeq: equality constraints matrix
% Aeq=[bus P; flows P; bus Q; mixed flows; offline units P; & Q]

n=length(Y);
% number of buses
wPg=size(PgU,2); dPg=size(PgU,3);
wPl=size(Pl,2); dPl=size(Pl,3);
% size("..",1)= bus, size("..",2)= bus components,
% size("..",3)= time interval
s_flows=length(nonzeros(triu(Y,1)));
% number of lines

if (dPg~=dPl)
    disp('horizons not matching');
    return;
end

lp_vars=1+(n-1)*(wPg+wPl)+n+2*s_flows;
% linearized OPF vars (Ps, Pl, Pg, deltas & active flows)
lp_eq=n+2*s_flows;
% linearized OPF equality constraints

lb=[]; ub=[]; ct=[];
for (i=1:1:dPg)
    lb=[lb;0];
    %lb=[lb;-1000];
    ub=[ub;1000];
    ct=[ct;cg(1,1,i)];
    % slack bus active power limits (mind lower as "0" or "-1000")
    % assumed as a horizontally slack bus (single generator - NO loads)
    lbp=[]; lbq=[]; ubp=[]; ubq=[]; ctp=[];
    for (j=2:1:n)
        lbp=[lbp;PgD(j,:,i)';Pl(j,:,i)'];
        ubp=[ubp;PgU(j,:,i)';Pl(j,:,i)'];
        lbq=[lbq;QgD(j,:,i)';(tan(acosf(j,:)).*Pl(j,:,i))'];
        ubq=[ubq;QgU(j,:,i)';(tan(acosf(j,:)).*Pl(j,:,i))'];
        ctp=[ctp;cg(j,:,i)';cl(j,:,i)'];
    end
    lb=[lb;lbp;0;-pi*ones(n-1,1);-1000*ones(2*s_flows,1)];
    % "0" = slack bus angle
    lb=[lb;0;lbq;1;0.9*ones(n-1,1);-1000*ones(2*s_flows,1)];
    % "0" slack bus reactive power limit, "1" slack bus voltage
    ub=[ub;ubp;0;pi*ones(n-1,1);1000*ones(2*s_flows,1)];
    % "0" = slack bus angle
    ub=[ub;1000;ubq;1;1.1*ones(n-1,1);1000*ones(2*s_flows,1)];
    % "1000" slack bus reactive power limit, "1" slack bus voltage
    ct=[ct;ctp];
    ct=[ct;zeros(size(ub,1)-size(ct,1),1)];
end

Aeq=[];
counter=1;

lenP=1+(n-1)*(wPg+wPl);
for i=1:1:n
    for j=1:1:i-1
        if (imag(Y(j,i))~=0)
            % Active power flow equations...
            Aeq(n+counter,lenP+j)=-imag(Y(j,i))+real(Y(j,i));
            Aeq(n+counter,lenP+i)=imag(Y(j,i))-real(Y(j,i));
            Aeq(n+s_flows+counter,lenP+j)=imag(Y(j,i))-real(Y(j,i));
            Aeq(n+s_flows+counter,lenP+i)=-imag(Y(j,i))+real(Y(j,i));
            Aeq(n+counter,lenP+n+counter)=1;
            Aeq(n+s_flows+counter,lenP+n+s_flows+counter)=1;
            % Mixed Active-Reactive power flow equations...
            Aeq(lp_eq+n+counter,lp_vars+lenP+j)=-(imag(Y(j,i))^2+real(Y(j,i))^2);
            Aeq(lp_eq+n+counter,lp_vars+lenP+i)=(imag(Y(j,i))^2+real(Y(j,i))^2);
            Aeq(lp_eq+n+counter,lenP+n+counter)=-real(Y(j,i));
            Aeq(lp_eq+n+counter,lp_vars+lenP+n+counter)=imag(Y(j,i));
            Aeq(lp_eq+n+s_flows+counter,lp_vars+lenP+i)=-(imag(Y(j,i))^2+real(Y(j,i))^2);
            Aeq(lp_eq+n+s_flows+counter,lp_vars+lenP+j)=(imag(Y(j,i))^2+real(Y(j,i))^2);
            Aeq(lp_eq+n+s_flows+counter,lenP+n+s_flows+counter)=-real(Y(j,i));
            Aeq(lp_eq+n+s_flows+counter,lp_vars+lenP+n+s_flows+counter)=imag(Y(j,i));
            
            % Active power flows summed to buses...
            % Reactive power flows summed to buses...
            Aeq(j,lenP+n+counter)=-1;
            Aeq(i,lenP+n+s_flows+counter)=-1;
            Aeq(lp_eq+j,lp_vars+lenP+n+counter)=-1;
            Aeq(lp_eq+i,lp_vars+lenP+n+s_flows+counter)=-1;
            if (i>1)
                Aeq(i,2+(i-2)*(wPg+wPl):1+(i-1)*(wPg+wPl))=[ones(1,wPg),-ones(1,wPl)];
                Aeq(lp_eq+i,lp_vars+2+(i-2)*(wPg+wPl):lp_vars+1+(i-1)*(wPg+wPl))=[ones(1,wPg),-ones(1,wPl)];
            end
            if (j>1)
                Aeq(lp_eq+j,lp_vars+2+(j-2)*(wPg+wPl):lp_vars+1+(j-1)*(wPg+wPl))=[ones(1,wPg),-ones(1,wPl)];
                Aeq(j,2+(j-2)*(wPg+wPl):1+(j-1)*(wPg+wPl))=[ones(1,wPg),-ones(1,wPl)];
            end
            counter=counter+1;
        end
    end
end
Aeq(1,1)=1; Aeq(lp_eq+1,lp_vars+1)=1;
% assumed as a horizontally slack bus (single generator - NO loads)
Aeqbuf=Aeq; bfr=1;
while bfr<dPg
    Aeq=sparse(blkdiag(Aeq,Aeqbuf));
    bfr=bfr+1;
end
%for (i=2:1:dPg)
%end

beq=zeros(size(Aeq,1),1);
    
% s=strcat('Aeqlt',datestr(now,'mmmmddyyyyHHMMSS'),'.txt');
% dlmwrite(s,horzcat([Aeq;lb';ub';ct']),'-append','delimiter','\t','newline','pc','precision','%.7f');

end