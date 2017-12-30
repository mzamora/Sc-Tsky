function make_streamer_input(date_i)
% Create input files for running Streamer using complete sounding profiles 
% for Sc clouds, and estimating cloud properties from a well-mixed assumption
% get Streamer at http://stratus.ssec.wisc.edu/streamer/
% 2017 Monica Zamora, SRAF at UCSD solar.ucsd.edu

%% Get sounding data & cloud properties
addpath('../Sc-utils/Soundings') % use the Sc-utils code
dmax=eomday(year(date_i),month(date_i)); %last day of month
pathName=['../Sc-utils/Soundings/raw/72293_',datestr(date_i,'yyyy_mm_'),'0100_',num2str(dmax),'12.csv']; %set source of data
ListofVar={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};
[p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(date_i,date_i+0.5,pathName,ListofVar);
p=double(p); z=double(z);
if isempty(z)
    return 
end
f=w>50; p(f)=[]; z(f)=[]; tp(f)=[]; RH(f)=[]; w(f)=[];
[~,~,hght_top,zi,eta_top,ni]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); %Xiaohui's code for inversion height
[~,~,hght_top2,~,eta_top2,~]=TMP_Inversion_Strength_Cal(-w,z/1000,z(1)); %Xiaohui's code for inversion height

if isempty(z)
    return 
end
if isnan(zi)
    fprintf('No inversion found \n')
    return
end
if zi==z(1)/1e3
    fprintf('Surface inversion \n')
    %return
end
if hght_top2<hght_top
    hght_top=hght_top2;
    eta_top=eta_top2;
end

%% Find cloud base from sounding
nocloud=0; %we assume initially that there'll be a cloud
try
    cloudpoints=find(RH>95);
    n_clouds=1;
    nb=cloudpoints(1);
    zb=z(cloudpoints(1));
    for ic=1:length(cloudpoints)-1
        if cloudpoints(ic+1)~=(cloudpoints(ic)+1)
            n_clouds=n_clouds+1;
            zb=z(cloudpoints(ic+1));
            nb=cloudpoints(ic+1);
        end
    end    
catch
    fprintf('No cloud layer \n')
    nocloud=1;
    zb=zi;
end

h=zi*1e3-zb; %cloud thickness in m
if h==0 
    fprintf('No cloud layer \n')
    nocloud=1;
end
Tcld=tp(ni); %temp at top of the cloud
PW=-trapz(p,w)/1000/9.81*10; %precipitable water in cm

%% Re-do vertical profile: clean as much as possible, but don't delete zi
i=1; notdone=1;
while notdone
    if (i+2)>=length(z)
        notdone=0;
        continue
    end
    tpi1=interp1([z(i),z(i+2)],[tp(i),tp(i+2)],z(i+1));
    if abs(tp(i+1)-tpi1)<0.5 && (z(i+1)~=zi*1e3)
        RH(i+1)=[]; p(i+1)=[]; z(i+1)=[]; tp(i+1)=[]; w(i+1)=[];
    end
    i=i+1;
end
    
%% Estimate tau,LWP,LWC from well mixed assumption
ni=find(z==round(zi*1e3,1),1); %update zi
if ni>1
    %%
    z2=z(1):1:zi*1e3;
    theta_l=trapz(z(1:ni),theta(1:ni))/(z(ni)-z2(1))*ones(size(z2));
    q_t=trapz(z(1:ni),w(1:ni))/(z2(end)-z2(1))*ones(size(z2))*1e-3; %in kg/kg
    p2=interp1(z(1:ni),p(1:ni),z2,'linear','extrap')*100; %in Pa
    [T,q_l,q_sat]=get_T_ql_qs(z2,theta_l,q_t,p2);

    zb2=z2(find(q_l>0,1)); %well mixed cloud base
    rho=p2./287.3./T;
    LWP=trapz(z2,rho.*q_l)*1000; %LWP in g/m2
    LWC=q_l(end)*rho(end)*1e3; %LWC at top
    tau=3*LWP/2/1e6/1e-5;
else
    nocloud=1; LWP=nan;
end
    %%
if LWP==0
    fprintf('No cloud in well-mixed approximation \n')
    return
end

if nocloud
    znew=z; RHnew=RH; tpnew=tp; pnew=p;
else
    %% Put 15 extra points in the cloud and 15 above it
    n=length(z);
    zCLD=linspace(zb,zi*1e3+h,30);
    znew=[z(z<zb);zCLD';z(z>zi*1e3+h)];
    RHnew=interp1(z,RH,znew);
    tpnew=interp1(z,tp,znew);
    pnew=interp1(z,p,znew);
end

%% Don't exceed 85 points for extending to 100km
if length(znew)>85
    znew=znew(1:85); RHnew=RHnew(1:85); pnew=pnew(1:85); tpnew=tpnew(1:85);
end

%% Create input file
filename_in=['NKX/sounding_' datestr(date_i,'yyyymmdd') '.inp'];
fileID = fopen(filename_in,'w');

%% Options section
fprintf(fileID,'OPTIONS \n');
fprintf(fileID,'.TRUE. \n'); %FLUXES compute fluxes?
fprintf(fileID,'.FALSE. \n'); %IR106 use the band 106(3.7um)?
fprintf(fileID,'.TRUE. \n'); %CLDFRC cloud forcing?
fprintf(fileID,'2 2 \n'); %NSTR* SW and LW number of streams
fprintf(fileID,'0 \n'); %NCOEF how many Legendre coeffs
fprintf(fileID,'.TRUE. \n'); %GASABS gaseous absorption?
fprintf(fileID,'.TRUE. \n'); %RAYLISHRT include Rayleigh scattering?
fprintf(fileID,'2 \n'); %ALBTYPE albedo control: model to be scaled by a number
fprintf(fileID,'4 \n'); %EMISSTYPE emissivity control: e for all bands
fprintf(fileID,'2 .TRUE. \n'); %Profile mid-lat summer, extend to 100km
fprintf(fileID,'4 1 \n'); %Aerosol model: maritime, background tropospheric aerosol
fprintf(fileID,'1 1 3 1 2 2 3\n'); %z(km) T(K) RH% O3 zi(km) h(km) λ(Streamer#)
fprintf(fileID,'4 \n'); %Output levels: all
fprintf(fileID,'.TRUE. \n'); %Descriptive output
fprintf(fileID,['../testio/NKX/sounding_',datestr(date_i,'yyyymmdd'),'.des \n']); %Output name
fprintf(fileID,'.TRUE. \n'); %YES user customized output
fprintf(fileID,['../testio/NKX/sounding_',datestr(date_i,'yyyymmdd'),'.out \n']); %
fprintf(fileID,' \n'); % No weights
fprintf(fileID,'.FALSE. \n'); % Won't read cloud optical properties
fprintf(fileID,' \n'); %
fprintf(fileID,' \n'); %

%% Case section
fprintf(fileID,'CASE \n'); %
fprintf(fileID,['Sounding NKX ',datestr(date_i,'yyyymmdd-hh'),'UTC \n']); %Title
fprintf(fileID,[datestr(date_i,'yy mm dd'),' 12.0 32.85 242.89 -99 \n']); %Date time and lat-long
fprintf(fileID,' \n'); % Nothing for FLUXES TRUE
fprintf(fileID,'1 129 \n'); %Bands
fprintf(fileID,' 0.05 1 1 1 \n'); %Albedo 0.05 w/1 sfc: water
fprintf(fileID,[num2str(tpnew(1)),' 1.0 \n']); %Tsrf and emissivity
if nocloud
    fprintf(fileID,'0 \n');
else
    fprintf(fileID,['1 1 0.9 ',num2str(Tcld),' ',num2str(zi),' ',num2str(h/1e3),' ',num2str(tau),' 1 0 10 ',num2str(LWC),' 999 999 999\n']); %Cloud info w/only 1 cloud properties
end
fprintf(fileID,' \n'); %No cloud overlap
fprintf(fileID,['1 1 1 1 0 0 ',num2str(length(znew)),' 1.0 1.0 1.0 1.0 1.0 1.0 \n']); % We read z,p,T,RH profiles
fprintf(fileID,'%.3f %.1f %.1f %.0f \n',[znew(end:-1:1)/1000,pnew(end:-1:1),tpnew(end:-1:1),RHnew(end:-1:1)]'); %put sounding profiles
fprintf(fileID,[num2str(znew(1)/1000),' ',num2str(PW),' \n']); %Sfc height and precipitable water(cm)

%% Close the file
fclose(fileID);

end
