function [ y ] = Sumzero( Nv,P)
Zi = [0.6099 0.0869 0.0691 0.0339 0.0378 0.0257 0.0212 0.0181 0.0601 0.0194 0.0121 0.0058];
MWi = [16.04 30.07 44.10 58.12 58.12 72.15 72.15 86.18 114.23 28.02 44.01 34.08 ];
Tbi = [94 303 416 471 491 542 557 610 743 109 194 331];		%Boiling temperature of components
Pci = [673 709 618 530 551 482 485 434 361 227 1073 672];	%Critical Pressure of components
Tci = [344 550 666 733 766 830 847 915 1024 492 548 1306];	%Critical Temperature of components
T = 150 +460;
Ppc = sum(Zi.*Pci); %critical pressure of sample
Tpc = sum(Zi.*Tci);	%critical temperature of sample
a = 1.2 + 4.5*(10^-4)*P + (1.5*(10^-9))*(P^2);
c = 0.89-1.7*(10^-4)*P - 3.5*(10^-8)*(P^2);
%Standings correlation for Ki
for i = 1:length(Zi)
    bi(i) = log10(Pci(i)/14.7)/((1/Tbi(i))-(1/Tci(i)));
    Fi(i) = bi(i)*((1/Tbi(i))-(1/T));
    ki(i) = (1/P)*(10^(a + c*Fi(i)));
    Y(i) = ((Zi(i)*(1-ki(i)))/((Nv*(ki(i)-1))+1));
end
y = sum(Y);

end

