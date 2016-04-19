from premia import pricings

results = set()

for m in pricings.all():
    for p in m.all():
        for meth in p.all():
            for x in meth.results():
               results.add(x)
            
for r in results:
   print r

