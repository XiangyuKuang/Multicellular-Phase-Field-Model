clear;
rehash toolboxcache

X=zeros(17,13)

for i=4:20
    %计算时间
    data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
    time=data.data;
    CD=data.textdata(2:end,1:2);
    [r1,col]=find(CD(:,2)=="ABa");
    time1=time(r1(:,1),1);
    [r2,col]=find(CD(:,2)=="ABp");
    time2=time(r2(:,1),1);
    [r3,col]=find(CD(:,2)=="ABal");
    time3=time(r3(:,1),1);
    [r4,col]=find(CD(:,2)=="ABar");
    time4=time(r4(:,1),1);
    [r5,col]=find(CD(:,2)=="ABpl");
    time5=time(r5(:,1),1);
    [r6,col]=find(CD(:,2)=="ABpr");
    time6=time(r6(:,1),1);
    [r7,col]=find(CD(:,2)=="EMS");
    time7=time(r7(:,1),1);
    [r8,col]=find(CD(:,2)=="P2");
    time8=time(r8(:,1),1);
    [r9,col]=find(CD(:,2)=="MS");
    time9=time(r9(:,1),1);
    [r10,col]=find(CD(:,2)=="E");
    time10=time(r10(:,1),1);
    [r11,col]=find(CD(:,2)=="C");
    time11=time(r11(:,1),1);
    [r12,col]=find(CD(:,2)=="P3");
    time12=time(r12(:,1),1);
    
    X(i-3,1)=max(time1)-min(time1)+1;
    X(i-3,2)=max(time2)-min(time2)+1;
    X(i-3,3)=max(time3)-min(time3)+1;
    X(i-3,4)=max(time4)-min(time4)+1;
    X(i-3,5)=max(time5)-min(time5)+1;
    X(i-3,6)=max(time6)-min(time6)+1;
    X(i-3,7)=max(time7)-min(time7)+1;
    X(i-3,8)=max(time8)-min(time8)+1;
    X(i-3,9)=max(time9)-min(time9)+1;
    X(i-3,10)=max(time10)-min(time10)+1;
    X(i-3,11)=max(time11)-min(time11)+1;
    X(i-3,12)=max(time12)-min(time12)+1;
    X(i-3,13)=max(time7);
    
end

for n=1:13
    a(n)=mean(X(:,n));
end



