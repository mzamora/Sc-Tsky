%read_all_NKX_outputs
addpath('../Sc-utils/Soundings') % use the Sc-utils code
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
            dates(ind)=date_i;
            qv_up3km(ind)=nan; qv_upper(ind)=nan;
            z_inv(ind)=nan; LWi_allz(ind)=nan; LWi_upperz(ind)=nan;
            
            try
            %% Get sounding data & cloud properties
            pathName=['../Sc-utils/Soundings/raw/72293_',datestr(date_i,'yyyy_mm_'),'0100_',num2str(dmax),'12.csv']; %set source of data
            ListofVar={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};
            [p,z,tp,~,RH,w,~,~,~,~,~]=Get_sounding_Var(date_i,date_i+0.5,pathName,ListofVar);
            
            [~,~,hght_top,zi,eta_top,ni]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); %Xiaohui's code for inversion height
            [~,~,hght_top2,~,eta_top2,~]=TMP_Inversion_Strength_Cal(-w,z/1000,z(1)); %Xiaohui's code for inversion height
            if hght_top2<hght_top
                hght_top=hght_top2;
                eta_top=eta_top2;
            end

            %% Compute upper values of moisture
            n_3km=find(z<3000,1,'last');
            qv_up3km(ind)=trapz(z(eta_top:n_3km),w(eta_top:n_3km))/(z(n_3km)-z(eta_top));
            qv_upper(ind)=trapz(z(eta_top:end),w(eta_top:end))/(z(end)-z(eta_top));
            catch
            end
            
            %% All profile
            filename=['NKX/sounding_' datestr(date_i,'yyyymmdd') '.out'];
            try
                [Z,P,T,WV,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
                z_inv(ind)=zi*1e3;
                LWi_allz(ind)=interp1(Z,LW_dn,zi);
            catch
            end
            %% Profile above inversion
            filename=['NKX2/sounding_' datestr(date_i,'yyyymmdd') '.out'];
            try
                [z2,p2,T2,wv2,O32,RH2,Aer2,LW_dn2,LW_up2]=read_streamer_output(filename);
                LWi_upperz(ind)=LW_dn2(end);
            catch
            end

            %% Plots
            try
            subplot(131); plot(LW_dn,Z); xlabel('LW ↓'); ylim([0 3]); ylabel('z'); hold on
            subplot(132); plot(LW_up,Z); xlabel('LW ↑'); ylim([0 3]); hold on
            title(datestr(date_i))
            subplot(133); plot(LW_up-LW_dn,Z); xlabel('Net LW'); ylim([0 3]); hold on
            legend('All z')
            catch
            end
            try
            subplot(131); plot(LW_dn2,z2); xlabel('LW ↓'); ylim([0 3]); ylabel('z')
            subplot(132); plot(LW_up2,z2); xlabel('LW ↑'); ylim([0 3])
            title(datestr(date_i))
            subplot(133); plot(LW_up2-LW_dn2,z2); xlabel('Net LW'); ylim([0 3])
            legend('z>z_i')
            catch
            end
            pause
            clf

            ind=ind+1;
        end
    end
end

%%
Tsky_allz=(LWi_allz/5.67e-8).^0.25;
Tsky_upperz=(LWi_upperz/5.67e-8).^0.25;

figure(1)
plot(LWi_allz,LWi_upperz,'.',[220 400],[220 400])
xlabel('LW_i full profile'); ylabel('LW_i upper profile')

figure(2)
plot(Tsky_allz,Tsky_upperz,'.',[250 295],[250 295])
xlabel('T_{sky} full profile'); ylabel('T_{sky} upper profile')

figure(3)
plot(Tsky_upperz,qv_up3km,'.')
xlabel('T_{sky}'); ylabel('q_v to 3km');

figure(4)
plot(Tsky_upperz,qv_upper,'.')
xlabel('T_{sky}'); ylabel('q_v upper');
%%
save('T_skies.mat','dates','qv_up3km','qv_upper','LWi_allz','LWi_upperz','Tsky_allz','Tsky_upperz')


%%
dates_far=dates(abs(LWi_allz-LWi_upperz)>20);
for i=1:length(dates_far)
    date_i=dates_far(i);
%%
    filename=['NKX/sounding_' datestr(date_i,'yyyymmdd') '.out'];
    [Z,P,T,WV,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
    
    filename=['NKX2/sounding_' datestr(date_i,'yyyymmdd') '.out'];
    [z2,p2,T2,wv2,O32,RH2,Aer2,LW_dn2,LW_up2]=read_streamer_output(filename);
    
    subplot(131); plot(LW_dn,Z,LW_dn2,z2); xlabel('LW ↓'); ylim([0 3]); ylabel('z')
    subplot(132); plot(LW_up,Z,LW_up2,z2); xlabel('LW ↑'); ylim([0 3])
    title(datestr(date_i))
    subplot(133); plot(LW_up-LW_dn,Z,LW_up2-LW_dn2,z2); xlabel('Net LW'); ylim([0 3])
    legend('All z','z>z_i')
    pause
end
