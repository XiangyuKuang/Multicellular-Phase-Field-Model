clear;
rehash toolboxcache
z=["ABa","ABp","P1"];
X=zeros(2,3,3)
for m=1:3
    for n=1:3
        for i=7
            %计算时间（ABa、ABp表面积同时存在的第一刻）
            t=8;

            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);


            [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
            if timeb(row1,t)
                X(1,m,n)=timeb(row1,t);
            end
        end
    end
end

for m=1:3
    for n=1:3
        for i=15
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row1,col1]=find(CD(:,2)=="ABa");
            time1=time(row1(:,1),1);
            [row2,col2]=find(CD(:,2)=="ABp");
            time2=time(row2(:,1),1);
            [row3,col2]=find(CD(:,2)=="P1");
            time3=time(row3(:,1),1);
            time1max=max(time1);
            time2max=max(time2);
            time3max=max(time3);
            ta=min(time1max,time2max);
            t=9;

            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
            if timeb(row1,t)
                X(2,m,n)=timeb(row1,t);
            end
        end
    end
end


answer=zeros(3,3)+0.0000000001;
answer1=zeros(3,3)+0.0000000001;
for m=1:3
    for n=1:3
        if sum(X(1:2,m,n))~=0;
            if find(X(1:2,m,n)==0)
                answer(m,n)=sum(X(1:2,m,n))./sum(sum(X(1:2,m,n)~=0));
                answer(n,m)=answer(m,n)
            else
                answer1(m,n)=sum(X(1:2,m,n))./sum(sum(X(1:2,m,n)~=0));
                answer1(n,m)=answer1(m,n)
            end
        end
    end
end

[r,c] = size(answer1);
x = 1:c;
y = 1:r;
[xx,yy] = meshgrid(x,y);
yy = flipud(yy);
% scatter(xx(:),yy(:),answer(:),"black");
hold on;
scatter(xx(:),yy(:),answer1(:),"black",'filled');
% 坐标轴美化
axis equal
% hTitle = title('Contact Area');
% set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'on', ...
    'TickDir', 'in', 'TickLength', [0 0], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
    'XTick', 0:1:c+1,...
    'XLim', [0 c+1],...
    'YTick', 0:1:r+1,...
    'YLim', [0 r+1],...
    'XTickLabel',{'','ABa','ABp','P1',''},...
    'YTickLabel',{'','P1','ABp','ABa',''})
set(gca,'FontSize',18,'Fontname', 'Arial');
a = 0:c+1
b = 0:c+1
b=-a+c+1;
plot(a,b,'k--','LineWidth',2.5)
