function [tp,q_l,q_sat]=get_T_ql_qs(z2,theta_l,q_t,p2)
% [tp,q_l,q_sat]=get_T_ql_qs(z,theta_l,q_t,p)
% This function creates T, q_sat and q_l profiles for well-mixes profiles
% Input: z[m], theta_l[K], q_t[kg/kg], p[Pa]
% Output: T[K], q_l[kg/kg]
% Monica Zamora, 2017. SRAF at UCSD solar.ucsd.edu

z2=z2-z2(1); %scale height to surface (if different to zero) 

%% Constants to be used
Rd=287.3; Lv=2.5e6; Cp=1015; rho=1.22;
exner=(p2/1e5).^(Rd/Cp);
a=610.94; b=17.652; c=243.04;

%% Get vertical profiles
for i=1:length(z2)
	% Properties assuming dry conditions
    t_dry(i)=theta_l(i)*exner(i);
    e_sat(i)=a*exp((t_dry(i)-273.15)*b/(t_dry(i)-273.15+c));
    q_sat(i)=0.622*e_sat(i)/(p2(i)-e_sat(i));
    tp(i)=t_dry(i);
    if q_sat(i)>q_t(i) %if q_sat>q_t, we are good: it was dry
        q_l(i)=0;
        continue
    else %if q_sat<q_t, we have liquid water
        eps=1;
        while eps>0.01 %we search iteratively for the real T of the moist mix
            t_new=theta_l(i)*exner(i)+Lv/Cp*(q_t(i)-q_sat(i));
            eps=abs(t_new-tp(i));
            tp(i)=0.5*t_new+0.5*tp(i); %update T as the average between old and new value to avoid quick jumps
            e_sat(i)=a*exp((tp(i)-273.15)*b/(tp(i)-273.15+c));
            q_sat(i)=0.622*e_sat(i)/(p2(i)-e_sat(i));
            q_l(i)=q_t(i)-q_sat(i);
        end
    end
end

end
