T = 150+460; 	%Temperature coversion
R = 10.73;		%Gas Constant in british units
Tsc = 520;		%Standard temperature
Psc = 14.7;		%Standard Pressure
Pin = 100:100:600; 
Go = 0.85;		%Specific gravity of stock-tank oil, water = 1
Gg = 0.65;		%Specific gravity of solution gas, air =1
Rs = 600;		%Gas solubility of oil, scf/STB
ncomp = 12; 	%Number of components

%Estimation of Nv: moles of liquid in vapor phase
	
for j = 1:length(Pin)
    Fun = @(Nv)Sumzero(Nv,Pin(j));
    nv(j) = fzero(Fun,0.5);
    nl(j) = 1-nv(j);
    kip(:,j) = Ki(Pin(j)); 	%Calculation of Ki at different pressures
end 


Zi = [0.6099 0.0869 0.0691 0.0339 0.0378 0.0257 0.0212 0.0181 0.0601 0.0194 0.0121 0.0058];	%Compressibility Factor of components
MWi = [16.04 30.07 44.10 58.12 58.12 72.15 72.15 86.18 114.23 28.02 44.01 34.08 ]; 			%Molecular Weights of components


k = 0;
j = 0;
%Calculation of x and y at different Pressures
for k = 1:length(Pin)
    for j = 1:ncomp
    y(j,k) = Zi(j)*kip(j,k)/(nl(k) + kip(j,k)*nv(k));
	%x(j,k) = kip(j,k)*y(j,k); Can calso be calculated by this step
    end
end
k = 0;
j = 0;
for k = 1:length(Pin)
    for j = 1:ncomp
		x(j,k) = Zi(j)/(nl(k) + kip(j,k)*nv(k));
   end
end
i = 0;

%Average Molecular weight calculation
for i = 1:length(Pin)
MWl(i) = sum(MWi.*x(:,i)');
end
i = 0;
for i = 1:length(Pin)
MWv(i) = sum(MWi.*y(:,i)');
end

i = 0;
%Average Z factor calculation
for i = 1:length(Pin)
z(i) = Zav(T,Pin(i));
end

%Calculation of density of vapour
for i = 1:length(Pin)
rhoV(i) = MWv(i)*Pin(i)/(z(i)*R*T);
end

%Calculation of density of liquid
rhoL = ((62.4*Go) + (0.0136*Rs*Gg))/(0.972 + 0.000147*((1.25*(T-460))+Rs*((Gg/Go)^0.5))^1.175);

%Calculation of volume of gas and liquid and GOR
for i = 1:length(Pin)
	Vvsc(i) = z(i)*nv(i)*R*Tsc/Psc;	%volume of vapour in standard conditions, scf
	Vl(i) = nl(i)*MWl(i)/rhoL;		%Volume of liquid phase, bbl
	GOR(i) =Vvsc(i)/Vl(i); 			%Gas-oil ratio scf/bbl
end

for i = 1:length(Pin)
API(i) = 141.5/(131.5+((16.0185*rhoV(i))/1000)); %API gravity in liquid phase
end

%Optional graph
%{ 
figure
subplot(2,1,1);
plot(GOR,Pin);
title('GOR vs Pressure');
xlabel('GOR');
ylabel('Pressure(Psia)');
subplot(2,1,2);
plot(Pin, API);
xlabel('API');
ylabel('Pressure(Psia)');
%[X,Y] = meshgrid(x(:,1),y(:,1));

%}

%Plots for the rest can be obtained by changing the index of x ie, 2 = for ethane and so on. 
title('Compositions of Methane')
plot(x(1,:), Pin);
hold on
plot(y(1,:),Pin);
hold on
xlabel('x,y')
ylabel('Pressure') 
plot(nv,Pin)
xlabel('Number of moles of vapour')
ylabel('Pressure')

hold on 
%plot(nl,Pin)