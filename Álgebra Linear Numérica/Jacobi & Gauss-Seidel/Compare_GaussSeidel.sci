//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//time1: tempo de compilação em GaussSeidel1 (com função inv)
//time2: tempo de compilação em GaussSeidel2 (com solução do sistema)
//////////////////////////////////////////////////////////////////////////
function [time1, time2]=Compare_GaussSeidel(n, E, M, type_norm)
    // matriz aleatoria
    A=n*eye(n,n)+rand(n,n);
    
    // vetores aleatorios
    b=rand(n,1);
    x0=rand(n,1);
    
    // calculando o tempo em GaussSeidel1
    tic();
    GaussSeidel1(A, b, x0, E, M, type_norm);
    time1 = toc();
    
    // calculando o tempo em GaussSeidel2
    tic();
    GaussSeidel2(A, b, x0, E, M, type_norm);
    time2 = toc();
endfunction
