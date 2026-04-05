function H=hank(x,k2,xpi,ypi,xi,yi,mi);
%function H=hank(x,k2,xpi,ypi,xi,yi,mi,ci);
%calcola la funzione di Hankel che costituisce la variabile 
%di integrazione per il calcolo dei coefficienti
%H=(besselh(0,2,(k2*sqrt([(xpi-x).^2]+([ypi-(yi+(mi*(x-xi)))].^2))))*ci);
H=besselj(0,k2*sqrt([(xpi-x).^2]+([ypi-(yi+(mi*(x-xi)))].^2)))-j*...
  bessely(0,k2*sqrt([(xpi-x).^2]+([ypi-(yi+(mi*(x-xi)))].^2)));