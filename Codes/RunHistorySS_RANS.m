clear;
clc;
close all;

nInterval = 10;
nRuns = 1;              % Enter the number of sim runs to index last run
nStart = 1;

axesNumfontsize = 15;

%% Data Read

allFolder = dir(append(pwd,'\postProcessing-1'));
allTimes = dir(append(pwd,'\postProcessing-1','\',allFolder(3).name));

for run=nStart:nRuns
    for fold=3:length(allFolder)
        if fold==3
            vars=string(allFolder(fold).name);
        else
            vars=[vars,string(allFolder(fold).name)];
        end

        for time=3:length(allTimes)
            if fold==3
                if time==3
                    times=string(allTimes(time).name);
                else
                    times=[times,string(allTimes(time).name)];
                end
            end
            file=append(pwd,'\postProcessing-',num2str(run),'\',allFolder(fold).name,'\',allTimes(time).name);
            allFiles=dir(file);
            for nfiles=3:length(allFiles)
                fileread=append(file,'\',allFiles(nfiles).name);

                %Find number of header lines in file
                fileID = fopen(fileread,'r');
                header="true";
                numHead=0;
                while header=="true"
                    currentLine = fgetl(fileID);
                    %Check if header line matches with the current line and store in iMatchOne:
                    if strfind(currentLine,"#")
                        numHead = numHead+1;
                    else
                        header="false";
                    end
                end

                ds{fold-2,run,nfiles-2} = datastore(fileread,'ReadVariableNames', true, 'NumHeaderLines', numHead);
            end
        end
    end
end
%% Plotting

% Plot mass flow
dsNum=[3, 4, 6, 9, 7, 11, 13];
%linetype=["-";".."];
plotColor=["red";"blue";"black";"green";"yellow"; "cyan";"magenta"];

figure
hold on
index=1;
for i=dsNum
    k=1;
    lastIteration=0;
    for t=nStart:nRuns
        data=readall(ds{i,t});
        count=nInterval;
        while count<length(data{:,1})
            x(k)=data{count,1}+lastIteration;
            y1(k)=sum(data{count-(nInterval-1):count,2})/nInterval;
            count=count+nInterval;
            k=k+1;
        end
        lastIteration=lastIteration+data{end,1};
    end

    plot(x,y1,'color',plotColor(index),'LineWidth',2);
    massflowIter=x;
    massflow(:,index)=y1;
    index=index+1;
end

xlabel("Iterations");
ylabel("Massflow (kg/s)");
grid on
grid minor
leg=legend(vars(dsNum),'Location','northeast','Box','off');
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');

% Plot net massflow
figure
hold on
netMassFlow=100.*(massflow(:,1)+massflow(:,2)+massflow(:,3)+massflow(:,4))./(massflow(:,4));
plot(massflowIter,netMassFlow,'LineWidth',2);
xlabel("Iterations");
ylabel("Massflow Error (%)");
grid on
grid minor
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');

%Plot change in massflow over iterations
%iter=floor(400/nInterval);
for i=1:length(massflow(:,1))
    if i==1
        massflowDelta(i,1:length(dsNum))=0;
    else
        massflowDelta(i,:)=abs(massflow(i,:)-massflow(i-1,:))./abs(massflow(i,:));
    end
end
plotColor=["red";"blue";"black";"green"];
figure
semilogy(massflowIter,massflowDelta(:,1),'color',plotColor(1),'LineWidth',2);
hold on
semilogy(massflowIter,massflowDelta(:,2),'color',plotColor(2),'LineWidth',2);
semilogy(massflowIter,massflowDelta(:,3),'color',plotColor(3),'LineWidth',2);
semilogy(massflowIter,massflowDelta(:,4),'color',plotColor(4),'LineWidth',2);
xlabel("Iterations");
ylabel("Massflow Error (%)");
grid on
grid minor
leg=legend(vars(dsNum),'Location','northeast','Box','off');
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');


%% Plot Max and min Temperature
dsNum=[1,2];
%linetype=["-";".."];
plotColor=["red";"blue"];
col=2; %temperature column
clear x;
clear y1;

figure
hold on
index=1;
for i=dsNum
    k=1;
    lastIteration=0;
    for t=nStart:nRuns
        data=readall(ds{i,t});
        count=nInterval;
        while count<length(data{:,1})
            x(k)=data{count,1}+lastIteration;
            y1(k)=data{count,col};
            count=count+nInterval;
            k=k+1;
        end
        lastIteration=lastIteration+data{end,1};
    end
    plot(x,y1,'color',plotColor(index),'LineWidth',2);
    index=index+1;
end

xlabel("Iterations");
ylabel("Temperatur (K)");
grid on
grid minor
leg=legend(vars(dsNum),'Location','northeast','Box','off');
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');

%% Plot Max and min Epsilon
dsNum=[1,2];
%linetype=["-";".."];
plotColor=["red";"blue"];
col=5; %epsilon column
clear x;
clear y1;

figure
hold on
index=1;
for i=dsNum
    k=1;
    lastIteration=0;
    for t=nStart:nRuns
        data=readall(ds{i,t});
        count=nInterval;
        while count<length(data{:,1})
            x(k)=data{count,1}+lastIteration;
            y1(k)=data{count,col};
            count=count+nInterval;
            k=k+1;
        end
        lastIteration=lastIteration+data{end,1};
    end
    plot(x,y1,'color',plotColor(index),'LineWidth',2);
    index=index+1;
end

xlabel("Iterations");
ylabel("Epsilon");
grid on
grid minor
leg=legend(vars(dsNum),'Location','northeast','Box','off');
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');



%% Plot Residuals

dsNum=[10];
%linetype=["-";".."];
plotColor=["#FF0000";"#0000FF";"#00FFFF";"#FF00FF";"#000000";"#FFFF00";"#00FF00";"#D95319"];
clear x;
clear y1;

figure
for i=dsNum
    k=1;
    lastIteration=0;
    for t=nStart:nRuns
        data=readall(ds{i,t});
        count=nInterval;
        while count<length(data{:,1})
            x(k)=data{count,1}+lastIteration;

            y1(k,1)=data{count,13};
            y1(k,2)=data{count,16};
            y1(k,3)=data{count,19};
            y1(k,4)=data{count,3};
            y1(k,5)=data{count,24};
            y1(k,6)=data{count,24};
            y1(k,7)=data{count,27};
            y1(k,8)=data{count,33};
            count=count+nInterval;
            k=k+1;

        end
        lastIteration=lastIteration+data{end,1};

        semilogy(x,y1(:,1),'color',plotColor(1),'LineWidth',2);
        hold on
        semilogy(x,y1(:,2),'color',plotColor(2),'LineWidth',2);
        semilogy(x,y1(:,3),'color',plotColor(3),'LineWidth',2);
        semilogy(x,y1(:,4),'color',plotColor(4),'LineWidth',2);
        semilogy(x,y1(:,5),'color',plotColor(5),'LineWidth',2);
        semilogy(x,y1(:,6),'color',plotColor(6),'LineWidth',2);
        semilogy(x,y1(:,7),'color',plotColor(7),'LineWidth',2);
        semilogy(x,y1(:,8),'color',plotColor(8),'LineWidth',2);
        %   leg=legend(["p","Ux","Uy", "h", "k", "epsilon","Rxx","Rxy", "Ryy"],'Location','northeast','Box','off');
        leg=legend(["Ux","Uy", "Uz","h", "epsilon","Rxx","Rxy", "Ryy"],'Location','northeast','Box','off');
    end
end

xlabel("Iterations");
ylabel("Residuals");
%axis([min(x) max(x) 1e-9 1e-2]);
grid on
grid minor
% leg=legend(["p","Ux","Uy","Uz","h","k","omega"],'Location','northeast','Box','off');
%leg=legend(["p","Ux","h","k","omega"],'Location','northeast','Box','off');
set(gca,'FontSize',axesNumfontsize, 'FontName', 'Times New Roman');
