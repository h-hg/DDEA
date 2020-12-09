
    
               fiteva=zeros(ps,1);%每一代的对应的个体被实际计算时为1，否则为0 
        minevafitpast=min(pastposfit);
    GPglobal;

    
        for mj=1:m
        k=ismember(p(mj,:),pastpos,'rows');
          if k==1
            pastid=find((ismember(p(mj,:),pastpos,'rows')));
            fiteva(mj,1)=1;%此处仅代表实际计算过，因为是赋值，所以并没有增加计算次数
            posfit(mj,1)=pastposfit(pastid(1),1);
            posmse(mj,1)=minmse;
          end
        end

       
        
        [FrontNo,MaxFNo] = NDSort(PopObj,m);%没有去掉期前的最优值，没关系，看fiteva如果已经为1则不计算
         difNo=unique(FrontNo);  %返回的是和FrontNo中一样的值，但是没有重复元素。产生的结果向量按升序排序。
         difNoL=length(difNo);%总共分了多少个面
         MinFNo=min(FrontNo);
         Fristpointsid=find (FrontNo==MinFNo);%找出在第一个面上点的下标
         Lastpointsid=find (FrontNo==MaxFNo);%找出在最后一个面上点的下标
        FristL=length(Fristpointsid);
        LastL=length(Lastpointsid);

            for ffid=1:FristL
              if fiteva(Fristpointsid(ffid),1)==0
              posfit(Fristpointsid(ffid),1)=evafit(p(Fristpointsid(ffid),:),problem);
                evacount=evacount+1;
                fiteva(Fristpointsid(ffid),1)=1;
                pastpos=[pastpos;p(Fristpointsid(ffid),:)];%%%更新历史数据
                pastposfit=[pastposfit;posfit(Fristpointsid(ffid),1)];
                posmse(Fristpointsid(ffid))=minmse;  
                minevafitnow = min(pastposfit);
               everyEVA=[everyEVA;evacount];
               everyGBEST=[everyGBEST;minevafitnow ];
 
             end
           end
      
             

for lfid=1:LastL
  if fiteva(Lastpointsid(lfid),1)==0

                posfit(Lastpointsid(lfid),1)=evafit(p(Lastpointsid(lfid),:),problem);
                evacount=evacount+1;
                fiteva(Lastpointsid(lfid),1)=1;
                pastpos=[pastpos;p(Lastpointsid(lfid),:)];%%%更新历史数据
                pastposfit=[pastposfit;posfit(Lastpointsid(lfid),1)];
               posmse(Lastpointsid(lfid))=minmse;  
                minevafitnow = min(pastposfit);
               everyEVA=[everyEVA;evacount];
               everyGBEST=[everyGBEST;minevafitnow ];
              
               
  end
end
           
        
       
      
      [posfit, gbestmid] = sort(posfit, 'descend');%%gbestmid为适应值由大到小排序，
      p=p(gbestmid,:);

      fiteva=fiteva(gbestmid,:);
    gbestval =minevafitnow;
    fprintf('Best fitness: %e\n',gbestval);
     fprintf('evacount: %e\n',evacount);
