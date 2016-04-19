from premia import pricings
from premia import models
from premia import options
from premia import assets

from time import time

for m in pricings.all():
    print "---processing model ", m
    for p in m.all():
        print "------processing pricing", p
        for meth in p.all():
            print "---------processing method", meth
            method = meth()
            model = method.model()()
            worst_time = 0
            for opt in meth.options():
                print "------------processing option", opt
                start = time()
                option = opt()
                try:
                    print option
                    print model
                    print method
                    print method(option, model)
                except Exception,exc:
                    print "EXCEPTION CAUGHT:", exc
                end = time()
                print end - start
                if worst_time < end - start:
                    worst_time = end - start
            print "Worst time:", worst_time
            
            if worst_time > 5: print "Worse5: ", method
            if worst_time > 10: print "Worse10: ", method            
            if worst_time > 20: print "Worse20: ", method            
            
