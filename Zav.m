function [ Zav ] = Zav(T,P)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
Zi = [0.6099 0.0869 0.0691 0.0339 0.0378 0.0257 0.0212 0.0181 0.0601 0.0194 0.0121 0.0058];
Pci = [673 709 618 530 551 482 485 434 361 227 1073 672];
Tci = [344 550 666 733 766 830 847 915 1024 492 548 1306];
Ppc = sum(Zi.*Pci);
Tpc = sum(Zi.*Tci);
Ppr = P/Ppc;
Tpr = T/Tpc;
A = 1.39*((Tpr-0.92)^0.5) - 0.36*Tpr -0.101;
B = (0.62-0.23*Tpr)*Ppr + ((0.066/(Tpr-0.86))-0.037)*Ppr^2 + (0.32*Ppr^6)/((10^(9*(Tpr-1))));
C = 0.132-0.32*log10(Tpr);
D = 10^(0.3106-0.49*Tpr+0.1824*Tpr^2);
Zav = A + ((1-A)/exp(B))+C*(Ppr^D);

end

