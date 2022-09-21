clear all
warning off
close;
Dim_set = [10, 20, 30, 50, 100];
%          F1   F2   F3   F4   F5   F6  F7
Fun_set = [101, 102, 103, 104, 105, 10, 19];
Reut_set = [];

global Dim Nd func_num

for i = 1:5
    Dim = Dim_set(1, i);
    for j = 1:7
        func_num = Fun_set(1, j);
        for k = 1:20
            Nd = k;
            output = BiS_SAHA;
            Reut(i, k) = output;
            str = ['Num_', num2str(i)];
            eval([str,'= Reut;']);
        end
    end
end