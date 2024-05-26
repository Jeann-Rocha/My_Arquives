//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz real n x n, simétrica
//Variáveis de saída:
//C: matriz n x 2, onde a i-ésima linha contém o centro e o raio do respectivo disco de Gerschgorin, nesta ordem
//////////////////////////////////////////////////////////////////////////

function [C] = Gerschgorin_real(A)
    [n]=size(A,1);
    C=zeros(n, 2);
    for (i=1:n)
        C(i,1)=A(i,i); //i-ésimo centro
        C(i,2)=sum(abs(A(i,:)))-abs(A(i,i)); //i-ésimo raio
    end
endfunction

//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz real n x n, simétrica
//x0: vetor, não nulo, a ser utilizado como aproximação inicial do autovetor dominante.
//epsilon: precisão a ser usada no critério de parada.
//alfa: valor do qual se deseja achar o autovalor de A mais próximo.
//M: número máximo de iterações.
//Variáveis de saída:
//lambdas: vetor com os autovalores encontrados
//////////////////////////////////////////////////////////////////////////

function [lambdas]=Find_eigenvalues(A,x0,epsilon,alfa,M)
    [n]=size(A,1);
    lambdas=zeros(n);
    for i=1:n
        [lambda1,x1,k,n_erro,time]=Potencia_deslocada_inversa(A,x0,epsilon,alfa(i),M); // autovalor mais próximo de alfa_i
        lambdas(i)=lambda1;
    end
endfunction
