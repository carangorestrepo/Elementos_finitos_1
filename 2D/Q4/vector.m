function A=vector(div,valores)
%div=0.1;
%valores=[0.5,1,1.2,1.7,3.2,3.5,4.4,5.5,6];
sizevalores=size(valores,2);
v1=valores(2:end)-valores(1:(sizevalores-1));
n=round(v1/div,0);
f=find(n==1);
n(f)=2;
divqa=v1./(n-1);
promedio=mean(divqa);

for i=1:(sizevalores-1)
    if divqa(i)>promedio
         n(i)=n(i)+1;
    end
end
n1=n-1;
%ns1=[1,cumsum(n)];
%ns=[1,cumsum(n1)];
%sumn=sum(n)-(sizevalores-1);
val=nan(sizevalores-1,max(n1));
for i=1:(sizevalores-1)
   valn=linspace(valores(i),valores(i+1),n(i));
   val(i,1:n1(i))=valn(1:n1(i));
end
n=(sizevalores-1)*max(n1);
B = reshape(val',[n,1]);
TF = isnan(B);
f=TF==0;
A=[B(f);valores(1,end)];
