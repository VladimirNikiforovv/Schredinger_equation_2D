clc
clearvars
cla
close all


L = 14;
T = 2.1;
dt = 0.001;
t = 0:dt:T;
dx = 0.27;
dy = 0.27;
x = -L/2:dx:L/2;
y = -L/2:dy:L/2;
Nx = length(x);
Ny = length(y);
TT = length(t);
a = dt/dx^2;
b = dt/dx^2;
N = Nx;


V = zeros(N,N);
for i = round(N/2)-round(N/16):round(N/2)
    %-round(N/8)
    for j = 1:round(N/2)-round(N/6)
        V(i,j) = 25;
    end
    for j = round(N/2)+round(N/6):N
        V(i,j) = 25;
    end
    for j = round(N/2)-round(N/16):N/2+round(N/16)
        V(i,j) = 25;
    end
end

mesh(x,y,V);

u = zeros(N^2,1);

for i = 1:N 
    for j = 0:N
        if j == 0
            u(i) = V(i,1);
        else
            u(i + Ny*j) = V(i,j);
        end       
    end
end


A = zeros(N^2,N^2);
for i = 1:N^2
    for j = 1:N^2
        if i == j
            A(i,j) = (-1i*a-1i*b-1i*dt*u(i)+1);
        elseif j == i+1  
            A(i,j) = 1i*a;
        elseif j == i-1   
            A(i,j) = 1i*a;
        elseif i + N == j
            A(i,j) = 1i*b;
        elseif i - N ==j
            A(i,j) = 1i*b;
        end
        if (mod((i),N) == 0 && mod((j-1),N) == 0) 
            A(i,j) = 0;
        end
        if (mod((i-1),N) == 0 && mod((j),N) == 0) 
            A(i,j) = 0;
        end
    end
end

% GausShl = zeros(N,N);
% for i = 1:Nx
%     for j = 1:Ny
%         GausShl(i,j) = exp(-((x(i) - 2 ).^2+(y(j) - 2).^2));
%     end
% end
% 
% mesh(GausShl);
% corelation -0.5.*(x(i) - 4 ).*(y(j) )
Psi = zeros(Nx,Ny,TT, 'like', complex(0,0));

for i = 1:Nx
    for j = 1:Ny
        Psi(i,j,1)= exp(-((x(i) - 3).^2+(y(j) - 0 ).^2)).*exp(1i.*(x(i)).*2.5 );
    end
end
% .*exp(1i.*(x(i)).*2.5 );


PsiVect = zeros(N^2,N^2);

for i = 1:Nx
    for j = 0:Ny-1
        if j == 0
            PsiVect(i,1) =  Psi(i,1,1);
        else
            PsiVect(i + Ny*j,1) =  Psi(i,j,1);
        end
    end
end




for g = 2:TT
    tic  
    PsiVect(:,g) = linsolve(A,PsiVect(:,g - 1));
    toc
    for i = 1:Nx
        for j = 0:N-1
            if j == 0
                Psi(i,1,g) = PsiVect(i,g);
            else
                Psi(i,j,g) = PsiVect(i + Ny*j,g);
            end            
        end
    end
end
% figure(2)
% mesh(x,y,abs(Psi(:,:,1).^2));
figure(2)

for k = 1:TT 
contour( abs(Psi(:,:,k)).^2,4);
pause(0.000001);
hold off
end

% for k = 500:2000
% figure(1)
% plot( abs(Psi(1,:,k)).^2)
% grid off
% pause(0.001);
% end

% dlmwrite ('C:/Users/nikiforov_va/Desktop/sar/Psi.txt', Psi, ' ')
% plot( abs(Psi(-1,:,1750)).^2)
% load('psi.mat')