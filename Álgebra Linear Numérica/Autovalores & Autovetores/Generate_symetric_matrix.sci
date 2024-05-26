//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//n: tamanho da matriz
//Variáveis de saída:
//A: matriz real n x n, simétrica
//////////////////////////////////////////////////////////////////////////

function [A] = Generate_symetric_matrix(n)
    A = zeros(n,n);
    for i=1:n
        for j=i:n
            A(i,j) = rand();
            A(j,i) = A(i, j); 
        end
    end
endfunction
