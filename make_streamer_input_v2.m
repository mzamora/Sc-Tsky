function make_streamer_input_v2(date_i)
% Create input files for running Streamer using sounding data for Sc clouds
% This version only gives out the properties above the cloud, since we are
% interested in the sky temperature
% get Streamer at http://stratus.ssec.wisc.edu/streamer/
% 2017 Monica Zamora, SRAF at UCSD solar.ucsd.edu

%% Get sounding data & cloud properties
addpath('../Sc-utils/Soundings') % use the Sc-utils code
dmax=eomday(year(date_i),month(date_i)); %last day of month
pathName=['../Sc-utils/Soundings/raw/72293_',datestr(date_i,'yyyy_mm_'),'0100_',num2str(dmax),'12.csv']; %set source of data
ListofVar={'PRES';'HGHT';'TEMP';'DWPT';'RELH';'MIXR';'WDIR';'WSPD';'THTA';'THTE';'THTV'};
[p z tp td RH w winddir windspeed theta theta_e theta_v]=Get_sounding_Var(date_i,date_i+0.5,pathName,ListofVar);
p=double(p); z=double(z);
f=w>50; p(f)=[]; z(f)=[]; tp(f)=[]; RH(f)=[]; w(f)=[];
[DT_max,numofInv,hght_top,zi,eta_top,ni]=TMP_Inversion_Strength_Cal(tp,z/1000,z(1)); %Xiaohui's code for inversion height

if isempty(z)
    fprintf('No data \n')
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
PW=-trapz(p,w)/1000/9.81*10;

if isempty(RH(RH>95 & z<3000))
    fprintf('No cloud in the BL \n')
    %return 
end

%% Clean arrays below the BL
z=z(ni:end); RH=RH(ni:end); p=p(ni:end); tp=tp(ni:end);

%% Re-do vertical profile: clean as much as possible
i=1; notdone=1;
while notdone
    if (i+2)>=length(z)
        notdone=0;
        continue
    end
    RHi1=interp1([z(i),z(i+2)],[RH(i),RH(i+2)],z(i+1)); %interpolation of neighbors
    if abs(RH(i+1)-RHi1)<0.5 %if the value in between and the interpolation are close
        RH(i+1)=[]; p(i+1)=[]; z(i+1)=[]; tp(i+1)=[]; %we take that point away
    end
    i=i+1;
end

%% Put 15 extra points above inversion (now surface) for better resolution
zup=linspace(z(1),zi*1.5e3,15)';

znew=[zup;z(z>1.5e3*zi)]; RHnew=interp1(z,RH,znew);
tpnew=interp1(z,tp,znew); pnew=interp1(z,p,znew);

%% limit the columns to 85 points (so it can extend the profiles to 100km)
if length(znew)>85
    znew=znew(1:85); RHnew=RHnew(1:85); pnew=pnew(1:85); tpnew=tpnew(1:85);
end
%% Create input file for streamer
filename_in=['NKX2/sounding_' datestr(date_i,'yyyymmdd') '.inp'];
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
fprintf(fileID,['../testio/NKX2/sounding_',datestr(date_i,'yyyymmdd'),'.des \n']); %Output name
fprintf(fileID,'.TRUE. \n'); %YES user customized output
fprintf(fileID,['../testio/NKX2/sounding_',datestr(date_i,'yyyymmdd'),'.out \n']); %
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
fprintf(fileID,'0.05 1 1 1 \n'); %Albedo 0.05 w/1 sfc: water
fprintf(fileID,[num2str(tpnew(1)),' 1.0 \n']); %Tsrf and emissivity
fprintf(fileID,'0 \n'); %No cloud this time
fprintf(fileID,' \n'); %No cloud overlap
fprintf(fileID,['1 1 1 1 0 0 ',num2str(length(znew)),' 1.0 1.0 1.0 1.0 1.0 1.0 \n']); % We read z,p,T,RH profiles
fprintf(fileID,'%.3f %.1f %.1f %.0f \n',[znew(end:-1:1)/1000,pnew(end:-1:1),tpnew(end:-1:1),RHnew(end:-1:1)]'); %put sounding profiles
fprintf(fileID,[num2str(znew(1)/1000),' ',num2str(PW),' \n']); %Sfc height and precipitable water(cm)

%% Close the file
fclose(fileID);

end