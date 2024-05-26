//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//R: Raio espectral da matriz do métodu de Gauss-Seidel.
//////////////////////////////////////////////////////////////////////////
function [R]=Raio_GaussSeidel(A)
    U=triu(A,1); // triangular superior de A
    invL=inv(A-U); // inversa de L
    
    R=max(abs(spec(-invL*(U)))); // raio espectral da matriz do métodu de Gauss-Seidel
endfunction
