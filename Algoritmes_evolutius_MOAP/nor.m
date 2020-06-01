
function ev = nor( x )
%Implimentacion para X inputs validos, solo si no hay ningun uno la salida
%sera 1
if sum(x) == 0
    ev = 1;
else
    ev = 0;
end

end

