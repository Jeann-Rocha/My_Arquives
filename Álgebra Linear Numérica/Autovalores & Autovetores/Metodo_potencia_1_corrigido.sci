//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz real n x n, diagonalizável, com autovalor dominante (lambda);
//x0: vetor, não nulo, a ser utilizado como aproximação inicial do autovetor dominante.
//epsilon: precisão a ser usada no critério de parada.
//M: número máximo de iterações.
//Variáveis de saída:
//lambda: autovalor dominante de A;
//x1: autovetor unitário (norma infinito) correspondente a lambda;
//k: número de iterações necessárias para a convergência;
//n_erro: norma infinito do erro
//Critério de parada: sendo erro=x1 – x0 (diferença entre dois iterados consecutivos), parar quando n_erro < épsilon ou k>M.
//////////////////////////////////////////////////////////////////////////

function [lambda,x1,k,n_erro,time] = Metodo_potencia_1_corrigido(A,x0,epsilon,M)
    tic();
    k=0;
    x0=x0/max(abs(x0));
    x1=A*x0; //aproximação do autovetor dominante
    n_erro=epsilon+1; //obriga a entrar no loop
    while (k<=M && n_erro>=epsilon)
        pos_max = find(abs(x1) == max(abs(x1)))(1); //posição da coordenada de maior módulo
        //aproximação do autovalor dominante com o sinal correto
        if (x1(pos_max) < 0)
            x1=-x1;
        end
        lambda=x1(pos_max);
        x1=x1/lambda;
        n_erro=norm((x1-x0),%inf);
        x0=x1;
        x1=A*x0;
        k=k+1;
    end
    time=toc();
endfunction
