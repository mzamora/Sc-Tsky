% create_all_NKX_inputs
% Creates inputs for streamer, from NKX soundings in the 2014-2017 summer months
% It will output them to a NKX folder
% Mónica Zamora, 2017. SRAF at UCSD solar.ucsd.edu

for yy=2014:2017
    for mm=5:9
        switch mm
            case {5,7,8}
                dmax=31;
            case {6,9}
                dmax=30;
        end
           
        for dd=1:dmax
            try
            make_streamer_input(datetime(yy,mm,dd,12,0,0));
            catch
            end
        end
    end
end
