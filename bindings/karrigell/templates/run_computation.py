import traceback

computation_timeout = 30

if g_errors <> []:
    run_computation = False

if 'Compute' not in REQUEST:
   run_computation = False

   
if run_computation and len(REQUEST) > n_elem(['m','f','o','meth',"Model_Size"]):
    try:
        def F(opt,model,q):
            begin = time()
            res = (method(opt,model))
            end = time()
            res.append(("Time", end - begin))
            q.put(res)

        def G(opt,model,q, itobj, itname, itlimit, itcount, itgetter, itsetter):
            initial = itgetter()            
            for i in range(itcount):
                begin = time()
                if itcount > 1:
                    x = initial + 1. * i / (itcount - 1) * (itlimit - initial)
                else:
                    x = initial
                    
                if type(initial) == type(1):
                    x = int(round(x))
                if type(initial) == type(1L):
                    x = long(round(x))
                
                itsetter(x)
                     
                res = (method(opt,model))
                end = time()
                res.insert(0, (itname, x))
                res.append(("Time", end - begin))
                q.put(res)
            
        res = []
        queue = Queue()
        if iterate_object == None:
            process = Process(target = F, args = (opt,model,queue))
        else:
            res_t = []
            process = Process(target = G, args = (opt,model,queue,iterate_object,iterate_name,iterate_to,int(_iteration_steps), iteration_getter, iteration_setter))
        process.start()
        try:
            if iterate_object == None:
                res.extend(queue.get(timeout=computation_timeout))
            else:
                begin = time()
                iterations = int(_iteration_steps)
                for i in range(iterations):
                    q = queue.get(timeout=begin + computation_timeout - time())
                    res.extend(q)
                    if res_t == []:
                        for k,v in q:
                            res_t.append((k,[v]))
                    else:
                        idx = 0
                        for k,v in q:
                            assert res_t[idx][0] == k
                            res_t[idx][1].append(v)
                            idx = idx + 1
                    
        except Empty, exc:
            if res == []: res = None
            if res_t == []: res_t = None
            process.terminate()
            raise Exception("Method has worked more than " + str(computation_timeout) + "s. Please try another parameter combination")
            
        process.join(timeout=computation_timeout)
    except Exception, exc:
        #print "exception caught:"
        #print traceback.format_exc()
        add_error(exc)

