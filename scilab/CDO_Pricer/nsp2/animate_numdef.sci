driver("X11")    
xselect()
xbasc()          
xset("pixmap",1)
realtimeinit(0.2);
realtime(0);
x=0.25:0.25:5;
rf=read('.losses_'+string(size(x,'*')/2),-1,2);
delta_r = rf(2,1)-rf(1,1);
xmax = rf(length(rf(:,1))-1,1);
ymax = max(rf(:,2));
if (tl == 1) rect = [0., 0., xmax, ymax];
else rect = [0., 0., xmax, ymax/delta_r]; end;

xset("wwpc")
for k=1:size(x,'*')
    s = string(x(k));
    if (length(s) == 1) s = s+'.00'; end;
    if (length(s) == 3) s = s+'0'; end;
    r=read('.losses_'+string(k-1),-1,2);
    realtime(k);
    xbasc();
    xset("wwpc")
    xtitle('t = '+s);
    if (tl == 1) plot2d3(r(:,1),r(:,2),2,rect=rect);
    else plot2d(r(:,1),r(:,2)/delta_r,2,rect=rect); end;
//    plot2d(r(:,1),r(:,2),2,rect=rect); 
    xset("wshow")
end;

