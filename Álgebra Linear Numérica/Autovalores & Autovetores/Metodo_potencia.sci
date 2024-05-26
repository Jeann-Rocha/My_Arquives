//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz real n x n, diagonalizável, com autovalor dominante (lambda);
//x0: vetor, não nulo, a ser utilizado como aproximação inicial do autovetor dominante.
//epsilon: precisão a ser usada no critério de parada.
//M: número máximo de iterações.
//Variáveis de saída:
//lambda: autovalor dominante de A;
//x1: autovetor unitário (norma infinito, no caso da Metodo_potencia_1 e, norma 2, no caso da Metodo_potencia_2) correspondente a lambda;
//k: número de iterações necessárias para a convergência;
//n_erro: norma infinito do erro
//Critério de parada: sendo erro=x1 – x0 (diferença entre dois iterados consecutivos), parar quando n_erro < épsilon ou k>M.
//////////////////////////////////////////////////////////////////////////

function [lambda,x1,k,n_erro,time] = Metodo_potencia_1(A,x0,epsilon,M)
    tic();
    k=0;
    x0=x0/max(abs(x0));
    x1=A*x0; //aproximação do autovetor dominante
    n_erro=epsilon+1; //obriga a entrar no loop
    while (k<=M && n_erro>=epsilon)
        lambda=max(abs(x1)); //aproximação do autovalor dominante
        x1=x1/lambda;
        n_erro=norm(x1-x0,%inf);
        x0=x1;
        x1=A*x0;
        k=k+1;
    end
    time=toc();
endfunction

function [lambda,x1,k,n_erro,time] = Metodo_potencia_2(A,x0,epsilon,M)
    tic();
    k=0;
    x0=x0/norm(x0,2);
    x1=A*x0; //aproximação do autovetor dominante
    n_erro=epsilon+1; //obriga a entrar no loop
    while (k<=M && n_erro>=epsilon)
        lambda=x1'*x0; //quociente de Rayleigh; x0 é unitário
        if (lambda < 0)
            x1=-x1; //mantém x1 com o mesmo sentido de x0
        end
        x1=x1/norm(x1,2);
        n_erro=norm(x1-x0,2);
        x0=x1;
        x1=A*x0;
        k=k+1;
    end
    time=toc();
endfunction
