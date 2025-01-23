function ARI=ari(c1,c2)
n=length(c1);
cp=crosstab(c1,c2);
cp(end+1,:)=sum(cp);
cp(:,end+1)=sum(cp,2);
%-----------------------------------------------
%º∆À„ARI
cp1=cp(1:end-1,1:end-1);
r0=sum(sum((cp1.*(cp1-1))./2));
cp2=cp(1:end-1,end)';
r1=sum((cp2.*(cp2-1))./2);
cp3=cp(end,1:end-1);
r2=sum((cp3.*(cp3-1))./2);
r3=(2*r1*r2)/(n*(n-1));
ARI=(r0-r3)/(0.5*(r1+r2)-r3);
end

