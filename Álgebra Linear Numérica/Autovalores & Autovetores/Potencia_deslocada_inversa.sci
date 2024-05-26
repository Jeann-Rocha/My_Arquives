//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz real n x n, diagonalizável;
//x0: vetor, não nulo, a ser utilizado como aproximação inicial do autovetor dominante.
//epsilon: precisão a ser usada no critério de parada.
//alfa: valor do qual se deseja achar o autovalor de A mais próximo.
//M: número máximo de iterações.
//Variáveis de saída:
//lambda: autovalor de A mais próximo de alfa;
//x1: autovetor unitário (norma 2) correspondente a lambda;
//k: número de iterações necessárias para a convergência;
//n_erro: norma infinito do erro
//Critério de parada: sendo erro=x1 – x0 (diferença entre dois iterados consecutivos), parar quando n_erro < épsilon ou k>M.
//////////////////////////////////////////////////////////////////////////

function [lambda1,x1,k,n_erro,time] = Potencia_deslocada_inversa(A,x0,epsilon,alfa,M)
    tic();
    k=0;
    [n]=size(A,1); //número de coordenadas
    I=eye(n,n); //matriz identidade n x n
    x0=x0/norm(x0,2);
    x1=x0;
    n_erro=epsilon+1; //obriga a entrar no loop
    while (k<=M && n_erro>=epsilon)
        [x1,C,P]=Gaussian_Elimination_4(A-alfa*I,x0); //resolve o sistema para achar x1
        lambda=x1'*x0; //quociente de Rayleigh; x1 é unitário
        x1=x1/norm(x1,2);
        if (lambda<0)
            x1=-x1; //mantém x1 com o mesmo sentido de x0
        end
        n_erro=norm(x1-x0,2);
        x0=x1;
        k=k+1;
    end
    lambda1=alfa+1/lambda;
    time=toc();
endfunction

//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//x: solução do sistema Ax=b (assumimos que tal solução existe).
//C: Seja A=LU a decomposição LU de A.
//Então C(i,j)=L(i,j) para i>j e C(i,j)=U(i,j) para j>=i.
//////////////////////////////////////////////////////////////////////////

function [x, C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    [n]=size(C,1);
    P=eye(n, n); //matriz de permutação inicializada (como identidade)
    for j=1:(n-1)
        //Caso o valor da posição pivô seja 0, a linha deve ser permutada com a primeira abaixo cujo elemento da coluna correspondente for diferente de 0
        [p_max, k] = max(abs(C(j:n, j))); //valor e posição do maior valor em módulo
        C([j,j+k-1],:) = C([j+k-1,j],:); //troca de linhas
        P([j,j+k-1],:) = P([j+k-1,j],:); //colocando elementos de troca de linhas na matriz
        
        //O pivô está na posição (j,j)
        for i=(j+1):n
        //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
        C(i,j)=C(i,j)/C(j,j);
        //Linha i  Linha i - C(i,j)*Linha j
        //Somente os elementos da diagonal ou acima dela são computados
        //(aqueles que compõem a matrix U)
        C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction
