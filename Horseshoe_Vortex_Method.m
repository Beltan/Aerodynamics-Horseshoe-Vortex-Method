function [cl,CL,CMle,CDi,y_cent,y,xf,xd] = HSVM(A,lambda,sweep,alpha,B,N,Uinf,graficos)

%- cl = distribución del coeficiente de sustentación
%- CL = coeficiente de sustentación del ala
%- CMle = coeficiente de momentos respecto el borde de ataque
%- CDi = coeficiente de resistencia inducida
%- y_cent = vector que define las posiciones del centro de cada panel
%- graficos (true/false) es para el caso que se desee visualizar el plot
%- A = Aspect Ratio o alargamiento
%- lambda = Taper Ratio o estrechamiento
%- sweep = ángulo de flecha (en rad)
%- alpha = ángulo de ataque (en rad)
%- B = envergadura del ala (en m)
%- N = numero de paneles que se desea tratar (N+1 son los puntos)
%- Uinf = velocidad de corriente libre (en m/s)

%Discretización de los bordes de los paneles
b=B/2; 
y=zeros(1,N+1);
xf=zeros(1,N+1);
xd=zeros(1,N+1);
c=zeros(1,N+1);
for i=(N/2+1):N+1
        y(i)=(i-(N/2+1))*2*b/N;
        xf(i)=tan(sweep)*y(i);
        xd(i)=((4/A)*((lambda-1)/(lambda+1))+tan(sweep))*y(i)+4*b/(A*(lambda+1));
        c(i)=xd(i)-xf(i);
end
for i=1:(N/2+1)
        y(i)=-y(N+2-i);
        xf(i)=xf(N+2-i);
        xd(i)=xd(N+2-i);
        c(i)=c(N+2-i);
end
%Discretización de los paneles
y_cent=zeros(1,N);
xf_cent=zeros(1,N);
xd_cent=zeros(1,N);
c_cent=zeros(1,N);
xc34=zeros(1,N);
xB=zeros(1,N);
xC=zeros(1,N);
xA=zeros(1,N);
xD=zeros(1,N);
for i=N/2:N
    y_cent(i)=(y(i+1)+y(i))/2;
    xf_cent(i)=tan(sweep)*y(i);
    xd_cent(i)=((4/A)*((lambda-1)/(lambda+1))+tan(sweep))*y(i)+4*b/(A*(lambda+1));
    c_cent(i)=xd_cent(i)-xf_cent(i);
    xc34(i)=c_cent(i)*3/4+xf_cent(i);
    xB(i)=xf(i)+c(i)/4;
    xC(i)=xf(i+1)+c(i+1)/4;
    xA(i)=xB(i)+40*b;
    xD(i)=xC(i)+40*b;
end
for i=1:N/2
    y_cent(i)=-y_cent(N+1-i);
    xf_cent(i)=xf_cent(N+1-i);
    xd_cent(i)=xd_cent(N+1-i);
    c_cent(i)=c_cent(N+1-i);
    xc34(i)=xc34(N+1-i);
    xB(i)=xB(N+1-i);
    xC(i)=xC(N+1-i);
    xA(i)=xA(N+1-i);
    xD(i)=xD(N+1-i);
end

%Cálculo de las velocidades inucidas, la matriz A, el RHS y la circulación
Vind=zeros(N,3,N);
A_matrix=zeros(N,N);
RHS=zeros(1,N);
wake=zeros(N,3,N);
for i=1:N
        for j=1:N
            Vind(i,:,j)=[0,0,0];
            Vind(i,:,j)=Vind(i,:,j)+induccion(xA(j),y(j),0,xB(j),y(j),0,xc34(i),y_cent(i),0);
            Vind(i,:,j)=Vind(i,:,j)+induccion(xC(j),y(j+1),0,xD(j),y(j+1),0,xc34(i),y_cent(i),0);
            wake(i,:,j)=Vind(i,:,j); %Necesareo para el CDi
            Vind(i,:,j)=Vind(i,:,j)+induccion(xB(j),y(j),0,xC(j),y(j+1),0,xc34(i),y_cent(i),0);
            A_matrix(i,j)=dot(Vind(i,:,j),[0,0,1]);
        end
        RHS(i)=-Uinf*sin(alpha);
end
gamma=(A_matrix)\(RHS.');

%Cálculo de la distribución de cl
S=(4*b^2)/A;
cl=zeros(1,N);
for i=1:N
    Si=(c(i+1)+c(i))*(y(i+1)-y(i))/2;
    cl(i)=(2*gamma(i)*(y(i+1)-y(i)))/(Uinf*Si);
end


%Cálculo del CL
sumgamma=0;
for i=1:N
    sumgamma=sumgamma+gamma(i)*(y(i+1)-y(i));
end
CL=(2/(Uinf*S))*sumgamma;

%Cálculp del CM en el borde de ataque
c_aero=0;
for i=N/2:N
    c_aero=c_aero+c_cent(i)^2*(y(i+1)-y(i));
end
c_aero=c_aero*2/S;
sumgamma=0;
for i=1:N
    sumgamma=sumgamma+gamma(i)*(xf_cent(i)+(c_cent(i)/4))*(y(i+1)-y(i))*cos(alpha);
end
CMle=-2*sumgamma/(Uinf*S*c_aero);

%Cálculo del CDi
sumgamma_alt=0;
for i=1:N
        sumgamma=0;
        for j=1:N
            sumgamma=sumgamma+gamma(j)*wake(i,3,j);
        end
        sumgamma_alt=sumgamma_alt+gamma(i)*(y(i+1)-y(i))*sumgamma;
end
CDi=-(2/((Uinf^2)*S))*sumgamma_alt;

%Graficación de los resultados obtenidos
if graficos
    figure(1)
    plot(y,xf,y,xd)
    if sweep >= 0
        axis([-b,b,-b/4,7*b/4])
    else
        axis([-b,b,-6*b/4,2*b/4])
    end  
    title('Forma del ala')
    figure(2)
    plot(y_cent,cl)
    title('Distribución del cl')
end