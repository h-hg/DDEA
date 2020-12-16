function [ obj,f,e,cost ] =MFB(POP,c,phi,problem_name)
% Usage: [ obj,f,e,cost ] =MFB(POP,c,phi,problem_name)
%
% Input:
% problem_name  - Multi-Fidelity Benchmark 
% c             - No. of Decision Variables
% POP           - Solution Set with c Decision Variables
% phi           - Fidelity Level [0,10000]
% 
%
% Output: 
% obj           - Evaluated Objective Values
% f             - Exact Objective Values
% e             - Errors
% cost          - Computational Cost

%
    %%%%    Authors:    Handing Wang, Yaochu Jin, John Doherty
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage/
    %%%%    DATE:       Oct 2017
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, John Doherty, A Generic Test Suite for
%Evolutionary Multi-Fidelity Optimization, IEEE Transactions on 
%Evolutionary Computation, Accepted, 10.1109/TEVC.2017.2758360

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
switch problem_name   
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
end


end