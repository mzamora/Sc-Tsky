%read_all_NKX_outputs
scutilsdir='../Sc-utils'; addpath([scutilsdir,'/Soundings']) % use the Sc-utils code
ind=1;

for yy=2014:2017
    for mm=5:9
        switch mm
            case {5,7,8}
                dmax=31;
            case {6,9}
                dmax=30;
        end
        for dd=1:dmax
            %% Whole profile
            date_i=datetime(yy,mm,dd,12,0,0)
            %%
            dates(ind)=date_i;
            qv_up3km(ind)=nan; qv_upper(ind)=nan;
            z_inv(ind)=nan; LWi_allz(ind)=nan; LWi_upperz(ind)=nan;
            
            try
            %% Get sounding data & inversion
            pathName=['../Sc-utils/Soundings/raw/72293_',datestr(date_i,'yyyy_mm_'),'0100_',num2str(dmax),'12.csv']; %set source of data
            ListofVar={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};
            [p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(date_i,date_i+0.5,pathName,ListofVar);
            p=double(p); z=double(z);
            f=w>50; p(f)=[]; z(f)=[]; tp(f)=[]; RH(f)=[]; w(f)=[];
            [~,~,zit,zi,nit,ni]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); %Xiaohui's code for inversion height
            [~,~,zit2,~,nit2,~]=TMP_Inversion_Strength_Cal(-w,z/1000,z(1)); %Xiaohui's code for inversion height
            if zit2<zit
                zit=zit2;
                nit=nit2;
            end

            %% Compute upper values of moisture
            n_3km=find(z<3000,1,'last');
            qv_up3km(ind)=trapz(z(nit:n_3km),w(nit:n_3km))/(z(n_3km)-z(nit));
            qv_upper(ind)=trapz(z(nit:end),w(nit:end))/(z(end)-z(nit));
            catch
                qv_up3km(ind)=nan;
                qv_upper(ind)=nan;
            end
            
            %% All profile
            filename=['NKX/sounding_' datestr(date_i,'yyyymmdd') '.out'];
            try
                [Z,P,T,WV,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
                z_inv(ind)=zi*1e3;
                LWi_allz(ind)=interp1(Z,LW_dn,zi,'nearest'); %,'extrap');
            catch
                try
                    [~,ni2]=min(abs(Z-zi));
                    LWi_allz(ind)=LW_dn(ni2);
                catch
                    LWi_allz(ind)=nan;
                end
            end
            %% Profile above inversion
            filename=['NKX2/sounding_' datestr(date_i,'yyyymmdd') '.out'];
            try
                [z2,p2,T2,wv2,O32,RH2,Aer2,LW_dn2,LW_up2]=read_streamer_output(filename);
                LWi_upperz(ind)=LW_dn2(end);
            catch
                LWi_upperz(ind)=nan;
            end

            %% Plots (just for checking)
%             try
%             subplot(131); plot(LW_dn,Z,'.-'); xlabel('LW ↓'); ylim([0 3]); ylabel('z'); hold on
%             text(150,zi*0.8,['LW_i all=',num2str(LWi_allz(ind))])
%             subplot(132); plot(LW_up,Z,'.-'); xlabel('LW ↑'); ylim([0 3]); hold on
%             title(datestr(date_i))
%             subplot(133); plot(LW_up-LW_dn,Z,'.-'); xlabel('Net LW'); ylim([0 3]); hold on
%             legend('All z')
%             catch
%             end
%             try
%             subplot(131); plot(LW_dn2,z2,'.-'); xlabel('LW ↓'); ylim([0 3]); ylabel('z')
%             text(150,zi*1.2,['LW_i upper=',num2str(LWi_upperz(ind))])
%             subplot(132); plot(LW_up2,z2,'.-'); xlabel('LW ↑'); ylim([0 3])
%             title(datestr(date_i))
%             subplot(133); plot(LW_up2-LW_dn2,z2,'.-'); xlabel('Net LW'); ylim([0 3])
%             legend('z>z_i')
%             catch
%             end
%             
%             pause
%             clf
            clear LW_dn LW_up LW_dn2 LW_up2 Z z2 zi
            ind=ind+1;
        end
    end
end

%% Get sky brightness temperature
Tsky_allz=(LWi_allz/5.67e-8).^0.25;
Tsky_upperz=(LWi_upperz/5.67e-8).^0.25;

%% Plots 1&2) Is the upper profile simulation giving the same result?
figure(1)
plot(LWi_allz,LWi_upperz,'.',[220 400],[220 400])
xlabel('LW_i full profile'); ylabel('LW_i upper profile')

figure(2)
plot(Tsky_allz,Tsky_upperz,'.',[250 295],[250 295])
xlabel('T_{sky} full profile'); ylabel('T_{sky} upper profile')

%% Plots 3&4) Do we see a trend between Tsky and qv metrics?
f=(~isnan(qv_up3km))&(~isnan(Tsky_upperz));
p_qv3km_Tsky=polyfit(qv_up3km(f),Tsky_upperz(f),1);
RSME_q3T=sqrt( sum((polyval(p_qv3km_Tsky,qv_up3km(f))-Tsky_upperz(f)).^2)/sum(f));
%q_qv3km_Tsky=polyfit(Tsky_upperz(f),qv_up3km(f),1); 
%RSME_Tq3=sqrt( sum((polyval(q_qv3km_Tsky,Tsky_upperz(f))-qv_up3km(f)).^2)/sum(f));

f=(~isnan(qv_up3km))&(~isnan(LWi_upperz));
p_qv3km_LWi=polyfit(qv_up3km(f),LWi_upperz(f),1);
RSME_q3L=sqrt( sum((polyval(p_qv3km_LWi,qv_up3km(f))-LWi_upperz(f)).^2)/sum(f));
%q_qv3km_LWi=polyfit(LWi_upperz(f),qv_up3km(f),1);
%RSME_Lq3=sqrt( sum((polyval(q_qv3km_LWi,LWi_upperz(f))-qv_up3km(f)).^2)/sum(f));

figure(3)
subplot(121)
plot(qv_up3km,LWi_upperz,'.',[0:0.1:12],polyval(p_qv3km_LWi,[0:0.1:12])) %,polyval(q_qv3km_LWi,[200:400]),[200:400])
ylabel('LW_i [W/m²]'); xlabel('q_v to 3km [g/kg]'); text(5,240,['LW_i=',num2str(p_qv3km_LWi(1),'%.2f'),'q_v+',num2str(p_qv3km_LWi(2),'%.2f')]);
legend('Data',['LW_i(q_v) RSME=',num2str(RSME_q3T)]) %,['q_v(LW_i) RSME=',num2str(RSME_Tq3)])
subplot(122)
plot(qv_up3km,Tsky_upperz,'.',[0:0.1:12],polyval(p_qv3km_Tsky,[0:0.1:12])) %,polyval(q_qv3km_Tsky,[250:290]),[250:290])
legend('Data',['T_{sky}(q_v) RSME=',num2str(RSME_q3L)]) %,['q_v(T_{sky}) RSME=',num2str(RSME_Tq3)])
ylabel('T_{sky} [K]'); xlabel('q_v to 3km [g/kg]'); text(5,257,['T_{sky}=',num2str(p_qv3km_Tsky(1),'%.2f'),'q_v+',num2str(p_qv3km_Tsky(2),'%.2f')]);

%%
f=(~isnan(qv_upper))&(~isnan(Tsky_upperz));
p_qvup_Tsky=polyfit(qv_upper(f),Tsky_upperz(f),1);
RSME_quT=sqrt( sum((polyval(p_qvup_Tsky,qv_upper(f))-Tsky_upperz(f)).^2)/sum(f));
%q_qvup_Tsky=polyfit(Tsky_upperz(f),qv_upper(f),1);
%RSME_Tqu=sqrt( sum((polyval(q_qvup_Tsky,Tsky_upperz(f))-qv_upper(f)).^2)/sum(f));

f=(~isnan(qv_upper))&(~isnan(LWi_upperz));
p_qvup_LWi=polyfit(qv_upper(f),LWi_upperz(f),1);
RSME_quL=sqrt( sum((polyval(p_qvup_LWi,qv_upper(f))-LWi_upperz(f)).^2)/sum(f));
%q_qvup_LWi=polyfit(LWi_upperz(f),qv_upper(f),1);
%RSME_Lqu=sqrt( sum((polyval(q_qvup_LWi,LWi_upperz(f))-qv_upper(f)).^2)/sum(f));

figure(4)
subplot(121)
plot(qv_upper,LWi_upperz,'.',[0:0.1:3],polyval(p_qvup_LWi,[0:0.1:3])) %,polyval(q_qvup_LWi,[200:400]),[200:400])
ylabel('LW_i [W/m²]'); xlabel('q_v upper [g/kg]'); text(1,240,['LW_i=',num2str(p_qvup_LWi(1),'%.2f'),'q_v+',num2str(p_qvup_LWi(2),'%.2f')]);
legend('Data',['LW_i(q_v) RSME=',num2str(RSME_quT)]) %,['q_v(LW_i) RSME=',num2str(RSME_Tqu)])
subplot(122)
plot(qv_upper,Tsky_upperz,'.',[0:0.1:3],polyval(p_qvup_Tsky,[0:0.1:3])) %,polyval(q_qvup_Tsky,[250:290]),[250:290])
ylabel('T_{sky} [K]'); xlabel('q_v upper [g/kg]'); text(1,257,['T_{sky}=',num2str(p_qvup_Tsky(1),'%.2f'),'q_v+',num2str(p_qvup_Tsky(2),'%.2f')]);
legend('Data',['T_{sky}(q_v) RSME=',num2str(RSME_quL)]) %,['q_v(T_{sky}) RSME=',num2str(RSME_Tqu)])

%% Save everything
save('results/T_skies.mat','dates','qv_up3km','qv_upper','LWi_allz','LWi_upperz','Tsky_allz','Tsky_upperz')

% load('results/T_skies.mat')

% %% For reviewing points far away
% dates_far=dates(abs(LWi_allz-LWi_upperz)>5);
% for i=1:length(dates_far)
%     date_i=dates_far(i);
% %%
%     filename=['NKX/sounding_' datestr(date_i,'yyyymmdd') '.out'];
%     [Z,P,T,WV,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
%     
%     filename=['NKX2/sounding_' datestr(date_i,'yyyymmdd') '.out'];
%     [z2,p2,T2,wv2,O32,RH2,Aer2,LW_dn2,LW_up2]=read_streamer_output(filename);
%     
%     subplot(131); plot(LW_dn,Z,LW_dn2,z2); xlabel('LW ↓'); ylim([0 3]); ylabel('z')
%     subplot(132); plot(LW_up,Z,LW_up2,z2); xlabel('LW ↑'); ylim([0 3])
%     title(datestr(date_i))
%     subplot(133); plot(LW_up-LW_dn,Z,LW_up2-LW_dn2,z2); xlabel('Net LW'); ylim([0 3])
%     legend('All z','z>z_i')
%     pause
% end
