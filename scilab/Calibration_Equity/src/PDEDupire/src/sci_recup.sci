function [stk,txt,top]=sci_recup()
RHS=[]
for k=1:rhs
  RHS=[stk(top)(1),RHS]
  top=top-1
end
top=top+1
//  w(i,1) is the ith output argument type
//  w(i,2) is the ith output argument row dimension
//  w(i,2) is the ith output argument column dimension
w=['?','?','?'
'?','?','?'
'?','?','?'
'?','?','?'
'?','?','?']
if lhs==1 then
  stk=list('recup'+rhsargs(RHS),'0',w(1,1),w(1,2),w(1,3))
else
  stk=list()
  for k=1:lhs
    stk(k)=list('recup'+rhsargs(RHS),'-1',w(k,1),w(k,2),w(k,3))
  end
end
