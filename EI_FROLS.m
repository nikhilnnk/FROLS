%% FROLS ALGORITHM TO FIT AN NARX MODEL
% Non Linear system
%          y_k :-0.605*y(k-1)-0.163*y^2(k-1)+0.588*u(k-1)-0.240*u(k-2)+e(k)
% where      e ~ Gaussian white noise
% Inputs:   ny : number of output delay terms
%           nu : number of input delay terms
%           ne : number of error terms
%            l : max degree
%            u : unifromily distribute input btw [-1 1]
% Outputs    c : selected terms (specified vector format)
%          not : number of selected terms 
%          phi : parameter values
%
%
%% Preliminaries
clear all;
clc;
ny=2                                                                         % nummber of y terms
nu=2                                                                         % nummber of u terms
ne=0;                                                                        % nummber of e terms
n=ny+nu+ne;                                                                  % total terms
l=3;                                                                         % max degree
ly=200;                                                                      % length 
M= factorial(n+l)/(factorial(n)*factorial(l))                                % Total number of possible terms
mu=0;                                                                        % noise average
sigma=0.2;                                                                   % noise standard deviation
e=sigma*randn(200,1)+mu;                                                     % Gaussian white noise
u=-1+2*rand(200,1);                                                          % uniformly distributed input
y=[0.1;0.5];                                                                 % initial output

for k=3:ly
    y(k)=(-0.605*y(k-1))-(0.163*y(k-2)^2)+(0.588*u(k-1))-(0.240*u(k-2))+e(k);%generating outputs
end

%% Generating all possible term sets
% C-stores all possible terms sets. --each row indicate a term
% and column values represent the powers corresponfing delay terms as u(k-1),u(k-2),y(k-1),y(k-2)
% eg: [1 0 0 2] indicate term = u(k-1)^1*y(k-2)^2
b=2;
C=zeros(M,n);                                                                % dimension (35,4)
for p=1:l 
    if p==1                                                                  % for power=1
      for i_1=1:n
          C(b,i_1)=C(b,i_1)+1;  
          b=b+1;
      end
    end
    if p==2                                                                  % for power=2
        for i_1=1:n
            for i_2=i_1:n
                C(b,i_1)=C(b,i_1)+1;
                C(b,i_2)=C(b,i_2)+1;
                b=b+1;
            end
        end
    end
    
    if p==3                                                                  % for power=3
        for i_1=1:n
            for i_2=i_1:n
                for i_3=i_2:n
                    C(b,i_1)=C(b,i_1)+1;
                    C(b,i_2)=C(b,i_2)+1;
                    C(b,i_3)=C(b,i_3)+1;  
                     b=b+1;
                end
            end
        end
    end     
   
end

%% creating D matrix -- Holds the values of all possible terms for each output
nm=max(nu,ny);                                                              % picking max delay
D=ones(ly-nm,M);
%size(D)
for i=nm+1:ly                                                               % Iteraing y
    for j=1:M                                                               % Iterating C
       k=1;
       for l=1:nu
       D(i-nm,j)=(D(i-nm,j))*(u(i-l)^(C(j,k)));                             % combining delay terms of input u
       k=k+1;
       end
       for m=1:ny
       D(i-nm,j)=(D(i-nm,j))*(y(i-m)^(C(j,k)));                             % combining delay terms of output y
       k=k+1;
       end
    end
end
size(D);

%% FROLS to Select the terms
Y=y(nm+1:ly,:);                         
sig=Y'*Y;                                                                    % sigma
sg=zeros(M,1);                                                               % vector to hold selected g_m
serr=zeros(M,1);                                                             % vector to hold selected err
%sl=zeros(M,1);
q=zeros(ly-nm,M);                                                            % vector to hold evaluated orthogonal vectors

for j=1:M
    err=zeros(M,1);
    g=zeros(M,1);
    if j==1                                                                  % loop to find first prominant term
        for i=1:M
            q(:,i)=D(:,i);
            qq=q(:,i)'*q(:,i);
            g(i)=(Y'*q(:,i))/qq;
            err(i)=(g(i)^2*qq)/sig;
        end
    else
        for i=1:M                                                            % loop to find second and remaining prominant terms
            if ~ismember(i,sl)
                p=D(:,i);
                pd=zeros(ly-nm,1);
                
                for r=1:size(sq,2)                                           % evaluating the subtracting term in orthogonalisation                                                             
                    qr=sq(:,r);
                    pd=pd+((p'*qr)/(qr'*qr))*qr;
                end
         
                q(:,i)=p-pd;                                                 % orthogonal vecctor q_m
                qq=q(:,i)'*q(:,i);                  
                g(i)=(Y'*q(:,i))/qq;
                err(i)=(g(i)^2*qq)/sig;                                      % evaluating err_m
            end
        end
    end       
[ERR,l]=max(err);                                                            % picking the term with max err and its index
serr(j)=ERR;                                                                 % store selected err value
sl(j)=l;                                                                     % store selected index of term
sg(j)=g(l);                                                                  % store selected g_m 
sq(:,j)=q(:,l);                                                              % store selected orthogonal vector

ESR=1-sum(serr);                                                             % termination parameter
if ESR<=0.05                                                                 % check termination condition
    break                                                                    % end the search if condition meets
end
end
sl;                                                                          % selected term index
c=C(sl,:)                                                                    % selected terms

%% Calculating the coefficients(parameters) using LS
for i=1:ly-nm
    X(i,:)=D(i,sl);
end
X;
not=size(c,1)                                                                % number of selected terms   
phi=inv(X'*X)*X'*Y                                                           % parameter values
          



