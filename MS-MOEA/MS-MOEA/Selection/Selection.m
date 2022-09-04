% An Adaptive Model Switch-based Surrogate-Assisted Evolutionary Algorithm 
% for Noisy Expensive Multi-Objective Optimization
%------------------------------- Reference --------------------------------
% N. Zheng, H. Wang, and B. Yuan, An Adaptive Model Switch-based 
% Surrogate-Assisted Evolutionary Algorithm for Noisy Expensive 
% Multi-Objective Optimization, Complex & Intelligent Systems, accepted, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Nan Zheng.
% Email: nanszheng@163.com
%当噪声影响作用小的时候，选择I（ipsl）+以及不确定性来选解
%当噪声影响作用大的时候，选择支配关系+pd以及不确定性来选解
function [Choose_Pop,Offspring_lastc,Offspring_lastd,sample_count]=Selection(Offspring_Covergence,Offspring_Diverity,MSE_Diverity,kappa,Offspring_lastc,Offspring_lastd,model,choose_Nflag)
        
        %% 根据I（ipsl）来进行适应度排序选解
         Offspring_savec=Offspring_Covergence;
         Offspring_saved=Offspring_Diverity;
         sample_count=zeros(1,3);
            if size(Offspring_lastc,1)==0
               OffCho = EnvironmentalSelection(Offspring_Covergence,5,kappa);
               sample_count(1)=1;
            else
                offobjc=Offspring_Covergence.objs;
                offdeclc=Offspring_lastc.decs;
                offobjd=Offspring_Diverity.objs;
                offdecld=Offspring_lastd.decs;
               offobjlc=cal_objectvalue(offdeclc,size(Offspring_Covergence.objs,2),choose_Nflag,model);
               offobjld=cal_objectvalue(offdecld,size(Offspring_Diverity.objs,2),choose_Nflag,model);
               Zmin = min([offobjc;offobjlc],[],1);
               h = cal_Convergence(offobjc,offobjlc,Zmin);
               if h==1
                  OffCho = EnvironmentalSelection(Offspring_Covergence,5,kappa);
                  sample_count(1)=1;
               else
                  if PD(offobjld) < PD(offobjd)
                     [OffCho,~]=Selection_uncertainty(Offspring_Diverity,MSE_Diverity,5);
                     sample_count(3)=1;
                  else
                    
                    Offspring_obj=Offspring_Diverity.decs;         
                    index = Paretoset(Offspring_obj); 
                    Offspring_Diverity=Offspring_Diverity(index');
                    if size(Offspring_Diverity,2)>5
                    OffCho = Selection_diverity(Offspring_Diverity,5); 
                    else 
                    OffCho = Offspring_Diverity;
                    end
                    sample_count(2)=1;
                  end
               end
               
            end

           Offspring_lastc=Offspring_savec;
           Offspring_lastd=Offspring_saved;
           Choose_Pop=[OffCho.decs];

end
