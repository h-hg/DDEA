function SDHJayamain()
clc;
for func_num=[1,3:30]
    dimension=10;
    nfes=[0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]*dimension*10000;
    savetimes=[]; 
    savemean=[]; 
    savestd=[]; 
    savegens=[];
    runs=1;
    Maxruns=11;
    savegen=[];
    titlenameMS=['SDHJayaMS_P',num2str(dimension),'D','.xlsx'];
    while(runs<Maxruns)
        titlenameTimes=['SDHJayaMS_',num2str(dimension),'D','_',num2str(func_num),'F','_',num2str(runs),'R','.txt'];
        titlenameGens=['SDHJayaGen_',num2str(dimension),'D','_',num2str(func_num),'F','_',num2str(runs),'R','.txt'];
        popsize=50;
        mini=-100*ones(1,dimension);
        maxi=100*ones(1,dimension); 
        x=zeros(popsize,dimension);
        for i=1:dimension
            x(:,i)=mini(i)+(maxi(i)-mini(i))*rand(popsize,1); 
        end
        f=fitness(x,func_num,popsize);
        historicalPop=x;
        xdata=x;
        fdata=f;
        gens=popsize;
%         maxgens=10000*dimension;
        maxgens=2000;
        m=1;
        Decision=1;
        Score=[1,0];
        Award=[0.5,0.5];
        Punishment=[(popsize/(popsize+gens*(5*dimension-1)))+popsize/(popsize+2),((gens*(5*dimension-1))/(popsize+gens*(dimension*dimension-1)))+2/(popsize+2)];
       
        savegen(runs,gens)=min(f);
        while(gens<maxgens)
            Stag=0;
            MaxStag=10;
            MaxSstag=20;
            Decision=1;
            Gbemin=min(fdata);
            Gafmin=min(fdata);
            while(Stag<MaxStag)
                bemin=min(fdata);
                afmin=min(fdata);
                switch Decision
                    case 1
                        xnew=[];
                        fnew=[];
                        if rand()<(cos((gens/maxgens)*pi)+1)*0.5
                            historicalPop=x;
                        end
                        xnew=popupdate(x,f);
                        if rand()<(cos((gens/maxgens)*2*pi)+1)*0.5
                            xnew=Crossover(popsize,dimension,x,historicalPop);
                        end
                        xnew=range(xnew,mini,maxi);
                        fnew=fitness(xnew,func_num,popsize);
                        xdata=[xdata;xnew];
                        fdata=[fdata,fnew];
                        for i=1:popsize
                            if(fnew(i)<f(i))
                                x(i,:)=xnew(i,:);
                                f(i)=fnew(i);
                            end
                        end
                        gens=gens+popsize;
                        bemin=afmin;
                        afmin=min(f);
                        if afmin<bemin
                            Score(Decision)=Score(Decision)+Award(Decision);
                        else
                            Score(Decision)=Score(Decision)- Punishment(Decision);
                        end
                    case 2
                        UModelPoint=[];
                        UModelfit=[];
                        BModelPoint=[];
                        BModelfit=[];
                        Uxdata=[];
                        Ufdata=[];
                        Txdata=[];
                        Tfdata=[];
                        [Uxdata,Ufdata,Txdata,Tfdata]=UdataSamp(xdata,fdata,dimension); 
                        UModelPoint=UdataModel(Uxdata,Ufdata,Txdata,Tfdata,xdata,fdata,mini,maxi,popsize,dimension,historicalPop,gens,maxgens);
                        UModelfit=fitness(UModelPoint,func_num,1);
                        gens=gens+1;
                        Uxdata=[Uxdata;UModelPoint];
                        Ufdata=[Ufdata,UModelfit];
                        [BModelPoint,fd]=BdataModel(Uxdata,Ufdata,Txdata,Tfdata,xdata,fdata,mini,maxi,popsize,MaxSstag,historicalPop,dimension,gens,maxgens); 
                        xdata=[xdata;UModelPoint];
                        fdata=[fdata,UModelfit];
                        if fd==1
                            BModelfit=fitness(BModelPoint,func_num,1);
                            gens=gens+1;
                            Uxdata=[Uxdata;BModelPoint];
                            Ufdata=[Ufdata,BModelfit];
                            xdata=[xdata;BModelPoint];
                            fdata=[fdata,BModelfit];
                            if BModelfit<min(fdata)
                                find=1;
                            else
                                find=0;
                            end
                        end
                        bemin=afmin;
                        afmin=min(fdata);
                        if afmin<bemin
                            Score(Decision)=Score(Decision)+Award(Decision);
                            [fdm,dindex]=min(fdata);
                            [ffm,findex]=min(f);
                            f(findex)=fdata(dindex);
                            x(findex,:)=xdata(findex,:);
                        else
                            Score(Decision)=Score(Decision)- Punishment(Decision);
                        end
                        mixfx=[fdata',xdata];
                        mixfx=sortrows(mixfx,1);
                        xdatasize=size(mixfx,1);
                        x=mixfx(1:popsize,2:end);
                        f=mixfx(1:popsize,1)';
                        
                end
                savegen(runs,gens)=min(f);
                Punishment=[(popsize/(popsize+gens*(5*dimension-1)))+popsize/(popsize+2),((gens*(5*dimension-1))/(popsize+gens*(dimension*dimension-1)))+2/(popsize+2)];
                if Score(1)>Score(2)
                    Decision=1;
                else
                    Decision=2;
                end
                Gbemin=Gafmin;
                Gafmin=min(fdata);
                if Gafmin<Gbemin
                    Stag=0;
                else
                    Stag=Stag+1;
                end
%                 if gens>=nfes(m)
%                     savegens(runs,m)=min(f);
%                     if gens>=nfes(m)
%                         m=m+1;
%                     end
%                 end
            end
            mixfx=[fdata',xdata];
            mixfx=sortrows(mixfx,1);
            xdatasize=size(mixfx,1);
            x=mixfx(1:5,2:end);
            f=mixfx(1:5,1)';
            r=randperm(xdatasize);
            for j=1:popsize-5
                x=[x;mixfx(r(j+5),2:end)];
                f=[f,mixfx(r(j+5),1)];
            end
            while(gens<maxgens)
                if rand()<(cos((gens/maxgens)*pi)+1)*0.5
                    historicalPop=x;
                end
                xnew=[];
                fnew=[];
                xnew=popupdate(x,f);
                if rand()<(cos((gens/maxgens)*2*pi)+1)*0.5
                    xnew=Crossover(popsize,dimension,x,historicalPop);
                end
                xnew=range(xnew,mini,maxi);
                gens=gens+popsize;
                fnew=fitness(xnew,func_num,popsize);
                for i=1:popsize
                    if(fnew(i)<f(i))
                        x(i,:)=xnew(i,:);
                        f(i)=fnew(i);
                    end
                end
                savegen(runs,gens)=min(f);
%                 if gens>=nfes(m)
%                     savegens(runs,m)=min(f);
%                     if gens>=nfes(m)
%                         m=m+1;
%                     end
%                 end
%                 min(f)
            end
        end
%         savetimes(1,runs)=min(f);
%         save(titlenameTimes,'savetimes','-ascii');
%         save(titlenameGens,'savegens','-ascii');
        runs=runs+1;
    end
    xlswrite( titlenameMS,savegen,func_num);
%     savemean=mean(savetimes,2);
%     savestd=std(savetimes,0,2); 
%     savemd=[savemean,savestd];
%     save(titlenameMS,'savemd','-ascii');
end
end