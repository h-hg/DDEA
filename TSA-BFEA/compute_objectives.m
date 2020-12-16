function [ obj,f,e,cost ] =compute_objectives(POP,c,phi,problem_name)
%计算目标值
%入口参数：POP-种群矩阵,c-变量个数，problem_name-测试问题
%出口参数：obj-目标函数值
switch problem_name
    case 'BPO1'
        load('-mat','modelpara1');
        load('-mat','modelpara2');
        load('coeff1.mat');
        n=size(POP,1);
        bu=[57.83	20.18	38.84	53.38	23.15	106.01	98.55	18.73	46.88	368.13	96.00	96.00	88.78	84.03	38.84];
        bd=[53.53	15.93	29.76	46.00	17.97	28.55	43.04	11.02	16.20	75.35	26.01	26.30	62.25	48.62	29.76];
        if phi==10000
            obj=evaluate_objective(POP, c, net1,net2);
            f=obj;
            e=0;
        else
            f=evaluate_objective(POP, c, net1,net2);
            POP=POP(:,1:c).*repmat(bu-bd,n,1)+repmat(bd,n,1);
            obj=evaluate_objectiveregression1(POP, c, mat_coef);
            e=obj-f;
        end
        cost=ones(n,1)*phi;
        
    case 'BPO2'
        load('-mat','modelpara1');
        load('-mat','modelpara2');
        load('coeff2.mat');
        n=size(POP,1);
        bu=[57.83	20.18	38.84	53.38	23.15	106.01	98.55	18.73	46.88	368.13	96.00	96.00	88.78	84.03	38.84];
        bd=[53.53	15.93	29.76	46.00	17.97	28.55	43.04	11.02	16.20	75.35	26.01	26.30	62.25	48.62	29.76];
        if phi==10000
            obj=evaluate_objective(POP, c, net1,net2);
            f=obj;
            e=0;
        else
            f=evaluate_objective(POP, c, net1,net2);
            POP=POP(:,1:c).*repmat(bu-bd,n,1)+repmat(bd,n,1);
            obj=evaluate_objectiveregression2(POP, c, mat_coef);
            e=obj-f;
        end
        cost=ones(n,1)*phi;
        
    case 'MFB1'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB2'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        theta=exp(-0.00025*phi);
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB3'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        if phi<1000
            theta=1-0.0002*phi;
        elseif phi<2000
            theta=0.8;
        elseif phi<3000
            theta=1.2-0.0002*phi;
        elseif phi<4000
            theta=0.6;
        elseif phi<5000
            theta=1.4-0.0002*phi;
        elseif phi<6000
            theta=0.4;
        elseif phi<7000
            theta=1.6-0.0002*phi;
        elseif phi<8000
            theta=0.2;
        elseif phi<9000
            theta=1.8-0.0002*phi;
        else
            theta=0;
        end
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*(0.001*phi).^4;
    case 'MFB4'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*(0.001*phi).^4;
    case 'MFB5'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        theta=exp(-0.00025*phi);
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB6'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        if phi==1000
            theta=0.8;
        else
            theta=0;
        end
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*phi;
        
    case 'MFB7'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        theta=1-0.0001*phi;
        psi=1-abs(POP(:,1:c));
        a=theta*ones(n,c).*psi;
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB8'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB9'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=c*exp(-0.0005*phi)*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        obj=f+e;
        cost=ones(n,1)*(0.001*phi).^4;
    case 'MFB10'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*0.1;
        mu=sum(1-abs(POP(:,1:c)),2)*sigma;
        e=randn(n,1)*sigma+mu;
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB11'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=exp(-0.0005*phi)*0.1;
        mu=sum(1-abs(POP(:,1:c)),2)*sigma;
        e=randn(n,1)*sigma+mu;
        obj=f+e;
        cost=ones(n,1)*(0.001*phi).^4;
    case 'MFB12'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB13'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=exp(-0.001*(phi+100));
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB1-Ackley'
        n=size(POP,1);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        
        a=20;
        b=0.2;
        cc=2*pi;
        POP(:,1:c)=POP(:,1:c)*32.768;
        f=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB1-Rastrigin'
        n=size(POP,1);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        f=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
        
        obj=f+e*10;
        cost=ones(n,1)*phi;
    case 'MFB1-Rosenbrock'
        n=size(POP,1);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        
        POP(:,1:c)=POP(:,1:c)*2.048;
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        f=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
        
        obj=f+e*100;
        cost=ones(n,1)*phi;
    case 'MFB1-Griewank'
        n=size(POP,1);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        
        POP(:,1:c)=POP(:,1:c)*600;
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        f=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB1-Ellipsoid'
        n=size(POP,1);
        theta=1-0.0001*phi;
        a=theta*ones(n,c);
        w=10*pi*theta*ones(n,c);
        b=0.5*pi*theta*ones(n,c);
        e=sum(a.*cos(w.*POP(:,1:c)+b+pi),2);
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        f=sum(P.*(POP(:,1:c).^2),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB8-Ackley'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        
        a=20;
        b=0.2;
        cc=2*pi;
        POP(:,1:c)=POP(:,1:c)*32.768;
        f=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB8-Rastrigin'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        f=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
        
        obj=f+e*10;
        cost=ones(n,1)*phi;
    case 'MFB8-Rosenbrock'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        
        POP(:,1:c)=POP(:,1:c)*2.048;
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        f=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
        
        obj=f+e*100;
        cost=ones(n,1)*phi;
    case 'MFB8-Griewank'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        
        POP(:,1:c)=POP(:,1:c)*600;
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        f=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB8-Ellipsoid'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        sigma=(1-0.0001*phi)*c*0.1;
        mu=0;
        e=randn(n,1)*sigma+mu;
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        f=sum(P.*(POP(:,1:c).^2),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB12-Ackley'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        
        a=20;
        b=0.2;
        cc=2*pi;
        POP(:,1:c)=POP(:,1:c)*32.768;
        f=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB12-Rastrigin'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        f=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
        
        obj=f+e*10;
        cost=ones(n,1)*phi;
    case 'MFB12-Rosenbrock'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        
        POP(:,1:c)=POP(:,1:c)*2.048;
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        f=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
        
        obj=f+e*100;
        cost=ones(n,1)*phi;
    case 'MFB12-Griewank'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        
        POP(:,1:c)=POP(:,1:c)*600;
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        f=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
    case 'MFB12-Ellipsoid'
        n=size(POP,1);
        f=c+sum(POP(:,1:c).^2-cos(10*pi*POP(:,1:c)),2);
        p=0.1*(1-0.0001*phi);
        e=zeros(n,1);
        for i=1:n
            if rand<=p
                e(i)=10*c;
            end
        end
        
        POP(:,1:c)=POP(:,1:c)*5.12;
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        f=sum(P.*(POP(:,1:c).^2),2);
        
        obj=f+e;
        cost=ones(n,1)*phi;
end


end