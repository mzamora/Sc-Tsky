for zcase=0:1:2
    filename=['DYCOMS/basic/DYCOMS',num2str(zcase,'%02d'),'.out'];
    [z,p,T,wv,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
    
    subplot(131); plot(LW_dn,z,'.-'); hold on
    subplot(132); plot(LW_up,z,'.-'); hold on
    subplot(133); plot(LW_up-LW_dn,z,'.-'); hold on
end
legend('Δz=105m','Δz=30m','Δz adap.')
subplot(131); xlabel('LW down'); ylabel('z[km]'); ylim([0 1.5])
subplot(132); xlabel('LW up'); ylim([0 1.5])
subplot(133); xlabel('Net LW'); ylim([0 1.5])

%% Basic vs extend up to 100km
    filename=['DYCOMS/basic/DYCOMS02.out'];
    [z,p,T,wv,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
    subplot(131); plot(LW_dn,z,'.-'); hold on
    subplot(132); plot(LW_up,z,'.-'); hold on
    subplot(133); plot(LW_up-LW_dn,z,'.-'); hold on
    
    filename=['DYCOMS/extended_up/DYCOMS02.out'];
    [z,p,T,wv,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
    subplot(131); plot(LW_dn,z,'.-'); hold on
    subplot(132); plot(LW_up,z,'.-'); hold on
    subplot(133); plot(LW_up-LW_dn,z,'.-'); hold on

legend('Basic','Extended')
subplot(131); xlabel('LW down'); ylabel('z[km]'); ylim([0 1.5])
subplot(132); xlabel('LW up'); ylim([0 1.5])
subplot(133); xlabel('Net LW'); ylim([0 1.5])

%% nstreams & ncoeff (all extended)
filenames={'DYCOMS/extended_up/DYCOMS02.out','DYCOMS/nstreams/DYCOMS4s2c.out','DYCOMS/nstreams/DYCOMS4s4c.out'};
for i=1:3
    filename=filenames{i};
    [z,p,T,wv,O3,RH,Aer,LW_dn,LW_up]=read_streamer_output(filename);
    subplot(131); plot(LW_dn,z,'.-'); hold on
    subplot(132); plot(LW_up,z,'.-'); hold on
    subplot(133); plot(LW_up-LW_dn,z,'.-'); hold on
end

legend('2 streams','4 streams, 2 coeff.','4 streams, 4 coeff.')
subplot(131); xlabel('LW down'); ylabel('z[km]'); ylim([0 1.5])
subplot(132); xlabel('LW up'); ylim([0 1.5])
subplot(133); xlabel('Net LW'); ylim([0 1.5])


