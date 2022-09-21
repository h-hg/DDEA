%%%global best
for r=1:20
    G_Optimal11(r,:) = G_Optimal(1,r,:);
end
save 'd:\FI\results\G_Optimal11.txt' G_Optimal11 -ascii

for r=1:20
    G_Optimal22(r,:) = G_Optimal(2,r,:);
end
save 'd:\FI\results\G_Optimal22.txt' G_Optimal22 -ascii
for r=1:20
    G_Optimal33(r,:) = G_Optimal(3,r,:);
end
save 'd:\FI\results\G_Optimal33.txt' G_Optimal33 -ascii
for r=1:20
    G_Optimal44(r,:) = G_Optimal(4,r,:);
end
save 'd:\FI\results\G_Optimal44.txt' G_Optimal44 -ascii
for r=1:20
    G_Optimal55(r,:) = G_Optimal(5,r,:);
end
save 'd:\FI\results\G_Optimal55.txt' G_Optimal -ascii
for r=1:20
    G_Optimal66(r,:) = G_Optimal(6,r,:);
end
save 'd:\FI\results\G_Optimal66.txt' G_Optimal66 -ascii
for r=1:20
    G_Optimal77(r,:) = G_Optimal(7,r,:);
end
save 'd:\FI\results\G_Optimal77.txt' G_Optimal77 -ascii
for r=1:20
    G_Optimal88(r,:) = G_Optimal(8,r,:);
end
save 'd:\FI\results\G_Optimal88.txt' G_Optimal88 -ascii

%%%%%%% Best for each evaluation
for r=1:20
    Iter_OptimalResult11(r,:) = Iter_OptimalResult(1,r,:);
end
save 'd:\FI\results\F11.txt' Iter_OptimalResult11 -ascii

for r=1:20
    Iter_OptimalResult22(r,:) = Iter_OptimalResult(2,r,:);
end
save 'd:\FI\results\F22.txt' Iter_OptimalResult22 -ascii
for r=1:20
    Iter_OptimalResult33(r,:) = Iter_OptimalResult(3,r,:);
end
save 'd:\FI\results\F33.txt' Iter_OptimalResult33 -ascii
for r=1:20
    Iter_OptimalResult44(r,:) = Iter_OptimalResult(4,r,:);
end
save 'd:\FI\results\F44.txt' Iter_OptimalResult44 -ascii
for r=1:20
    Iter_OptimalResult55(r,:) = Iter_OptimalResult(5,r,:);
end
save 'd:\FI\results\F55.txt' Iter_OptimalResult55 -ascii
for r=1:20
    Iter_OptimalResult66(r,:) = Iter_OptimalResult(6,r,:);
end
save 'd:\FI\results\F66.txt' Iter_OptimalResult66 -ascii
for r=1:20
    Iter_OptimalResult77(r,:) = Iter_OptimalResult(7,r,:);
end
save 'd:\FI\results\F77.txt' Iter_OptimalResult77 -ascii
for r=1:20
    Iter_OptimalResult88(r,:) = Iter_OptimalResult(8,r,:);
end
save 'd:\FI\results\F88.txt' Iter_OptimalResult88 -ascii

%%%% Time
for r=1:20
    time11(r,:) = Time(1,r,:);
end
save 'd:\FI\results\time11.txt' time11 -ascii

for r=1:20
    time22(r,:) = Time(2,r,:);
end
save 'd:\FI\results\time22.txt' time22 -ascii
for r=1:20
    time33(r,:) = Time(3,r,:);
end
save 'd:\FI\results\time33.txt' time33 -ascii
for r=1:20
    time44(r,:) = Time(4,r,:);
end
save 'd:\FI\results\time44.txt' time44 -ascii
for r=1:20
    time55(r,:) = Time(5,r,:);
end
save 'd:\FI\results\time55.txt' time55 -ascii
for r=1:20
    time66(r,:) = Time(6,r,:);
end
save 'd:\FI\results\time66.txt' time66 -ascii
for r=1:20
    time77(r,:) = Time(7,r,:);
end
save 'd:\FI\results\time77.txt' time77 -ascii
for r=1:20
    time88(r,:) = Time(8,r,:);
end
save 'd:\FI\results\time88.txt' time88 -ascii