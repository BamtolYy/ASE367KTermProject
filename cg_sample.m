function [cg, mass_total, moment] = cg_sample()
%point mass mean locations
x0 =[0,0.3,0.7,1.1,1.2,2.2,2.3,2.9,3.1,3.3,3.4,4,4.2,4.3,4.9,6,6.5,6.7,6.9,7.1,7.6];
%uniform sample location errors
x_err = 0.002*(1-2*rand(size(x0)));
%point mass means
mass0 = [0,20,25,25,20,25,70,1,10,15,1,20,120,1,50,5,1,5,5,10,10];
%uniform sample mass errors
mass_err = 0.1*(1-2*rand(size(x0)));
%add errors to means
x = x0 + x_err;
mass = mass0 + mass_err;
%compute properties
mass_total = sum(mass);
cg = sum(mass.*x)/mass_total;
moment = sum( mass.*(x-cg).^2 );
end
