clear;
rehash toolboxcache
z=["ABal","ABar","ABpl","ABpr","MS","E","C","P3"];
X=zeros(17,8,8)
for m=1:8
    for n=1:8
        for i=4:20
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [r1,col]=find(CD(:,2)=="ABal");
            time1=time(r1(:,1),1);
            [r2,col]=find(CD(:,2)=="ABar");
            time2=time(r2(:,1),1);
            [r3,col]=find(CD(:,2)=="ABpl");
            time3=time(r3(:,1),1);
            [r4,col]=find(CD(:,2)=="ABpr");
            time4=time(r4(:,1),1);
            time1max=max(time1);
            time2max=max(time2);
            time3max=max(time3);
            time4max=max(time4);
            ta=min(time1max,time2max);
            tb=min(time3max,time4max);
            t=min(ta,tb);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
            dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
            timea=dataa.data;
            surface=dataa.textdata;
            [~,cola]=find(surface(1,:)=="ABal");
            [~,colb]=find(surface(1,:)=="ABar");
            [~,colc]=find(surface(1,:)=="ABpl");
            [row,cold]=find(surface(1,:)=="ABpr");
            s0=timea(t,cola(1,1));
            s1=timea(t,colb(1,1));
            s2=timea(t,colc(1,1));
            s3=timea(t,cold(1,1));
            S=[s0,s1,s2,s3];
            A=isnan(S);
            B=sum(A);
            if B==0
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                    X(i-3,m,n)=timeb(row1,t);
                end
            end
        end
    end
end

answer=zeros(8,8)+0.0000000001;
answer1=zeros(8,8)+0.000000001;
for m=1:8
    for n=1:8
        if sum(X(1:17,m,n))~=0;
            if find(X(1:17,m,n)==0)
                answer(m,n)=sum(X(1:17,m,n))./sum(sum(X(1:17,m,n)~=0));
                answer(n,m)=answer(m,n)
            else
                answer1(m,n)=sum(X(1:17,m,n))./sum(sum(X(1:17,m,n)~=0));
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
         'XTickLabel',{'','ABal','ABar','ABpl','ABpr','MS','E','C','P3',''},...
         'YTickLabel',{'','P3','C','E','MS','ABpr','ABpl','ABar','ABal',''})
set(gca,'FontSize',18,'Fontname', 'Arial');
a = 0:c+1
b = 0:c+1
b=-a+c+1;
plot(a,b,'k--','LineWidth',2.5)






