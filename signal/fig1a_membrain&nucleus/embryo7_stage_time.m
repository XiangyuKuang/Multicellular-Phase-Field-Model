clear;
a=[];
load('E:\data\cellname.mat');

for i=7
    data=importdata(['E:\data\CD\CD',num2str(i,'%02d'),'.csv'] ) ;
    time=data.data;
    CD=data.textdata(2:end,1:2);
    t=[];
    for j=1:7
        z=string(Cellname{1,j});
        initial_time=[];
        for m=1:length(z)
            [row,~]=find(CD(:,2)==z(m));
            time1=time(row(:,1),1);
            time1max=max(time1);
            initial_time=[initial_time,time1max];
            t0=min(initial_time);
        end
        t=[t,t0];
    end
end