
function [ev,w] = nueva_moap_0706(A,w)
% Cost function for the landscape example
% [elev] = testfunction(loc)
table = [0 0 0; 0 0 1; 0 1 0; 0 1 1 ; 1 0 0; 1 0 1; 1 1 0;1 1 1];
f1 = [0 0 0 1 0 1 0 1];
f2 = [1 1 0 0 0 0 0 1];
f3 = [0 1 1 0 0 0 0 1];
f4 = [0 1 1 1 0 0 0 1];
sz = size(A);
dh = 0; %distancia de hamming: diferencia entre componentes de un vector--> cuanta mas diferencia mayor coste genera
cxs = 0; 
ev = 0;
for op = 1:8%totes les opcions de convinatoria  
n = zeros(1,(sz(2)+3));
n(1:3) = table(op,:); %valor de los tres primeros nodos, que corresponden a los inputs
    for k = 1:sz(2) % k = 1,1r nodo
        input = [];
                for i = 1:sz(1) %miro quienes son los inputs
                         if A(i,k) == 1
                             input = [input,n(i)];
                             if op ==8
                                 cxs = cxs + 1;
                             end
                         end
                end
                n(k + 3) = nor(input);
    end
     if f3(op)~= n(sz(2)+3) 
        
        w(op) = w(op) + 0.0001; %cada iteracion va subiendo el peso
    else
        w(op) = w(op) - 0.0001;
        
     end

    if op == 8
      dh = dh + 1.5*sqrt((n(sz(2)+3)-f4(op))^2);
%     elseif op == 2
%     dh = dh + 1.5*sqrt((n(sz(2)+3)-f3(op))^2);  
%     elseif op == 3
%     dh = dh + 2.5*sqrt((n(sz(2)+3)-f3(op))^2);
%     elseif op == 4
%     dh = dh + 2.5*sqrt((n(sz(2)+3)-f3(op))^2);
    else
     dh = dh + sqrt((n(sz(2)+3)-f4(op))^2);
    end
 
end
 ev = 0.99*dh + 0.01*(cxs+ sz(2));
end







%     if f3(op)~= n(sz(2)+3) 
%         
%         w(op) = w(op) + 0.0001; %cada iteracion va subiendo el peso
%     else
%         w(op) = w(op) - 0.00001;
%         
%     end





     