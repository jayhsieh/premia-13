HOME=getenv("HOME") + "/";
PREMIADIR=HOME + "devel/premia/trunk-dev";


exec(PREMIADIR + '/nsp/libpremia/loader.sce');
premia_init()

pb_list = PBLIST;

t=cputime();
Lpb=pb_list;
for pb=Lpb'
    printf("job %s\n",pb)
    load(pb);
    P.compute[];
end
t=t - cputime();
printf ("cpu : %f\n", t);
