//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//R: matriz do métodu de Jacobi
//////////////////////////////////////////////////////////////////////////
function [R]=Raio_Jacobi(A)
    d=diag(A); // vetor diagonal de A
    D=diag(d); // matriz da diagonal de A
    invD=diag(1./d); // inversa de D
    
    R=max(abs(spec(-invD*(A-D)))); // raio espectral da matriz do métodu de Jacobi
endfunction
