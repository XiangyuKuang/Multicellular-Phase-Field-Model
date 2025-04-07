clear;
rehash toolboxcache
z=["AB","P1"];
X=zeros(2,2,2)
for m=1:2
    for n=1:2
        for i=7
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row1,col1]=find(CD(:,2)=="AB");
            time1=time(row1(:,1),1);
            [row2,col2]=find(CD(:,2)=="P1");
            time2=time(row2(:,1),1);
            time1max=max(time1);
            time2max=max(time2);
            t=min(time1max,time2max);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
%             dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
%             timea=dataa.data;
%             surface=dataa.textdata;
%             [~,cola]=find(surface(1,:)=="ABa");
%             [~,colb]=find(surface(1,:)=="ABp");
%             s0=timea(t,cola(1,1));
%             s1=timea(t,colb(1,1));
%             S=[s0,s1];
%             A=isnan(S);
%             B=sum(A);
%             if B==0
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                    X(1,m,n)=timeb(row1,t);
                end
%             end
        end
    end
end

for m=1:2
    for n=1:2
        for i=15
            %计算时间
            data=importdata(['E:\elegans\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
            time=data.data;
            CD=data.textdata(2:end,1:2);
            [row1,col1]=find(CD(:,2)=="AB");
            time1=time(row1(:,1),1);
            [row2,col2]=find(CD(:,2)=="P1");
            time2=time(row2(:,1),1);
            time1max=max(time1);
            time2max=max(time2);
            t=min(time1max,time2max);
            
            %导入接触面积
            datab=importdata(['E:\elegans\data\Stat\Sample',num2str(i,'%02d'),'_Stat.csv']) ;
            timeb=datab.data;
            Stat=datab.textdata(:,1:2);
            
            %导入表面积
%             dataa=importdata(['E:\elegans\data\surface/Sample',num2str(i,'%02d'),'_surface.csv']) ;
%             timea=dataa.data;
%             surface=dataa.textdata;
%             [~,cola]=find(surface(1,:)=="ABa");
%             [~,colb]=find(surface(1,:)=="ABp");
%             s0=timea(t,cola(1,1));
%             s1=timea(t,colb(1,1));
%             S=[s0,s1];
%             A=isnan(S);
%             B=sum(A);
%             if B==0
                [row1,co1]=find((Stat(:,1)==z(m))&(Stat(:,2)==z(n)));
                if timeb(row1,t)
                    X(2,m,n)=timeb(row1,t);
                end
%             end
        end
    end
end


answer=zeros(2,2)+0.0000000001;
answer1=zeros(2,2)+0.0000000001;
for m=1:2
    for n=1:2
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
         'XTickLabel',{'','AB','P1',''},...
         'YTickLabel',{'','P1','AB',''})
set(gca,'FontSize',18,'Fontname', 'Arial');
a = 0:c+1
b = 0:c+1
b=-a+c+1;
plot(a,b,'k--','LineWidth',2.5)
