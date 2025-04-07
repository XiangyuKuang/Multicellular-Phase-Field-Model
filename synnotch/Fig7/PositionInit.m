function rc = PositionInit(N,L,r)
rc = rand(N,3).*(L-2*r)+r;
k = 0.1;
t = 1;
contact = ones(N,1);
delta = 0.01;
while t<20000&&(any(contact)||any([rc-r<0,rc+r-L>0],'all'))
  rc_temp = rc;
  for n = 1:N
    d = vecnorm(rc_temp-rc_temp(n,:),2,2);
    ContIdx = d<2*r;
    ContIdx(n) = 0;
    x1 = L-rc_temp(n,:)-r-delta;
    x2 = rc_temp(n,:)-r-delta;
%     rc(n,:) = rc(n,:)-k*(sum(rc_temp(ContIdx,:)-rc_temp(n,:),1)+2*abs(x1).*(x1<0)-2*abs(x2).*(x2<0)); % npj F=k(ri-rj)
    rc(n,:) = rc(n,:)-k*(sum( (1-2*r./d(ContIdx)).*(rc_temp(n,:)-rc_temp(ContIdx,:)) ,1)+2*abs(x1).*(x1<0)-2*abs(x2).*(x2<0)); % F=k(r-r0)v
    contact(n) = any(ContIdx);
  end
  t = t+1;
end