z0=50;
z_L=25-j*20;
fc=6e9;

X_L=imag(z_L);
R_L=real(z_L);
if R_L >= z0
    B1=(X_L+sqrt(R_L/z0)*sqrt(R_L^2+X_L^2-z0*R_L))/(R_L^2+X_L^2);
    B2=(X_L-sqrt(R_L/z0)*sqrt(R_L^2+X_L^2-z0*R_L))/(R_L^2+X_L^2);
    X1=1/B1+X_L*z0/R_L-z0/(B1*R_L);
    X2=1/B2+X_L*z0/R_L-z0/(B2*R_L);
end

if R_L < z0
    B1=sqrt((z0-R_L)/R_L)/z0;
    B2=-sqrt((z0-R_L)/R_L)/z0;
    X1=sqrt(R_L*(z0-R_L))-X_L;
    X2=-sqrt(R_L*(z0-R_L))-X_L;
end

f=0:1e4:12e9;
if B1>=0
    B1prime=B1*f/fc;
end
if B1<0
    B1prime=B1*fc./f;
end

if X1>=0
    X1prime=X1*f/fc;
end
if X1<0
    X1prime=X1*fc./f;
end

if B2>=0
    B2prime=B2*f/fc;
end
if B2<0
    B2prime=B2*fc./f;
end

if X2>=0
    X2prime=X2*f/fc;
end
if X2<0
    X2prime=X2*fc./f;
end

if R_L >= z0
    z_Lprime1=j*X1prime+1./(j*B1prime+ones(1,length(f))/(R_L+j*X_L));
    z_Lprime2=j*X2prime+1./(j*B2prime+ones(1,length(f))/(R_L+j*X_L));
end

if R_L < z0
    z_Lprime1=1./(j*B1prime+ones(1,length(f))./(R_L+j*(X1prime+X_L)));
    z_Lprime2=1./(j*B2prime+ones(1,length(f))./(R_L+j*(X2prime+X_L)));
end
gamma1=(z_Lprime1-ones(1,length(f))*z0)./(z_Lprime1+ones(1,length(f))*z0);
gamma2=(z_Lprime2-ones(1,length(f))*z0)./(z_Lprime2+ones(1,length(f))*z0);
hold on;
plot(f/1e9,abs(gamma1), 'DisplayName', 'Solution 1');
plot(f/1e9,abs(gamma2), 'red','DisplayName', 'Solution 2');
xlabel('Frquency(GHz)');
ylabel('Gamma');
legend;
hold off;
