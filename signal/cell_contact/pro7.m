clear;
rehash toolboxcache
z=["ABal","ABar","ABpl","ABpr","MS","E","P2"];
X=zeros(15,7,7)
for m=1:7
    for n=1:7
        for i=4:12
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row,col]=find(CD(:,2)=="P2");
            time1=time(row(:,1),1);
            t=max(time1);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
%             dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
%             timea=dataa.data;
%             surface=dataa.textdata;
%             [row,cola]=find(surface(1,:)=="P3");
%             [row,colb]=find(surface(1,:)=="C");
%             s0=timea(t,cola(1,1));
%             s1=timea(t,colb(1,1));
%             if isnan(s0) && isnan(s1)
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                    X(i-3,m,n)=timeb(row1,t);
                end
%             end
        end
    end
end

for m=1:7
    for n=1:7
        for i=14:15
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row,col]=find(CD(:,2)=="P2");
            time1=time(row(:,1),1);
            t=max(time1);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
%             dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
%             timea=dataa.data;
%             surface=dataa.textdata;
%             [row,cola]=find(surface(1,:)=="P3");
%             [row,colb]=find(surface(1,:)=="C");
%             s0=timea(t,cola(1,1));
%             s1=timea(t,colb(1,1));
%             if isnan(s0) && isnan(s1)
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                 X(i-4,m,n)=timeb(row1,t);
                end
%             end
        end
    end
end

for m=1:7
    for n=1:7
        for i=17:20
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row,col]=find(CD(:,2)=="P2");
            time1=time(row(:,1),1);
            t=max(time1);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
%             dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
%             timea=dataa.data;
%             surface=dataa.textdata;
%             [row,cola]=find(surface(1,:)=="P3");
%             [row,colb]=find(surface(1,:)=="C");
%             s0=timea(t,cola(1,1));
%             s1=timea(t,colb(1,1));
%             if isnan(s0) && isnan(s1)
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                X(i-5,m,n)=timeb(row1,t);
                end
%             end
        end
    end
end

answer=zeros(7,7)+0.0000000001;
answer1=zeros(7,7)+0.000000001;
for m=1:7
    for n=1:7
        if sum(X(1:15,m,n))~=0;
            if find(X(1:15,m,n)==0)
                answer(m,n)=sum(X(1:15,m,n))./sum(sum(X(1:15,m,n)~=0));
                answer(n,m)=answer(m,n)
            else
                answer1(m,n)=sum(X(1:15,m,n))./sum(sum(X(1:15,m,n)~=0));
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
%scatter(xx(:),yy(:),answer(:),"black");
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
         'XTickLabel',{'','ABal','ABar','ABpl','ABpr','MS','E','P2',''},...
         'YTickLabel',{'','P2','E','MS','ABpr','ABpl','ABar','ABal',''})
set(gca,'FontSize',18,'Fontname', 'Arial');
a = 0:c+1
b = 0:c+1
b=-a+c+1;
plot(a,b,'k--','LineWidth',2.5)