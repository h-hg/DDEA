function model_ensemble = build_ensemble_model(DataBase)

global Dim

%% 1 
sub_DataBase = DataBase;  
S = sub_DataBase(:, 1:Dim); Y = sub_DataBase(:, Dim + 1);
flag = 'cubic';
[lambda, gamma]=RBF(S,Y,flag); trainx = S;
FUN = @(x) RBF_eval(x,trainx,lambda,gamma,flag);
model_ensemble{1} = FUN;

%% 2
sub_DataBase = DataBase;  
S = sub_DataBase(:, 1:Dim); Y = sub_DataBase(:, Dim + 1);
flag = 'invmultiquad';
[lambda, gamma]=RBF(S,Y,flag); trainx = S;
FUN = @(x) RBF_eval(x,trainx,lambda,gamma,flag);
model_ensemble{2} = FUN;





