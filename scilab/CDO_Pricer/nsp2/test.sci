ncomp=100;
nominal=0.01*ones(100,1);
dates=[0.25:0.25:5];
tranches=[0,0.03,0.06,0.1,1.];
intensity=nominal;
xrates=[0.,5.];
yrates=[0.03,0.15];
[price,pl,dl]=price_cdo(ncomp, nominal, dates, tranches, intensity, xrates, yrates, 1, 0.4, 1, 0.03, 2, [2, 100, 1000]);
