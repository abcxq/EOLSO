function [fitness,X_error,U_duan]=fobj_RC_R0(x,U,I,ocv,R0,model_RC)  %U是电压，I是电流，都是行向量


if  model_RC==1    %  一阶RC模型
    
%     R0= x(1);
    R1= x(1);
    C1= x(2);

    L=size(I,1);
    X_error=zeros(L,1);
    U_duan=zeros(L,1);
    Ts=1;
    u1=0;
    %  error=0;

     for k=1:L

         u1 = exp(-Ts/(R1*C1))*u1+R1*(1-exp(-Ts/(R1*C1)))*I(k);
         UL = ocv(k)+I(k)*R0+u1;
         X_error(k) = U(k)-UL;
         U_duan(k) = UL;
         
     end

elseif  model_RC==2     %  二阶RC模型
    
%     R0= x(1);
    R1= x(1);
    R2= x(2);
    C1= x(3);
    C2= x(4);

    L=size(I,1);
    X_error=zeros(L,1);
    U_duan=zeros(L,1);
    Ts=1;
    u1=0;
    u2=0;
    %  error=0;

     for k=1:L

         u1 = exp(-Ts/(R1*C1))*u1+R1*(1-exp(-Ts/(R1*C1)))*I(k);
         u2 = exp(-Ts/(R2*C2))*u2+R2*(1-exp(-Ts/(R2*C2)))*I(k);
         UL = ocv(k)+I(k)*R0+u1+u2;
         X_error(k) = U(k)-UL;
         U_duan(k) = UL;
         
     end
     
elseif  model_RC==3     %  三阶RC模型
    
%     R0= x(1);
    R1= x(1);
    R2= x(2);
    R3= x(3);
    C1= x(4);
    C2= x(5);
    C3= x(6);

    L=size(I,1);
    X_error=zeros(L,1);
    U_duan=zeros(L,1);
    Ts=1;
    u1=0;
    u2=0;
    u3=0;
    %  error=0;

     for k=1:L

         u1 = exp(-Ts/(R1*C1))*u1+R1*(1-exp(-Ts/(R1*C1)))*I(k);
         u2 = exp(-Ts/(R2*C2))*u2+R1*(1-exp(-Ts/(R2*C2)))*I(k);
         u3 = exp(-Ts/(R3*C3))*u3+R1*(1-exp(-Ts/(R3*C3)))*I(k);
         UL = ocv(k)+I(k)*R0+u1+u2+u3;
         X_error(k) = U(k)-UL;
         U_duan(k) = UL;
         
     end
     
end
     
%  以均方根误差为目标函数     
fitness= sqrt(mean(abs(X_error).*abs(X_error)));  %   均方根误差

 
 
 