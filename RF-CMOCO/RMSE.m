function [ E ] = RMSE( POP,TPOP )
E=mean(abs(POP-TPOP),1);
end

