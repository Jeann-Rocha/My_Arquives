//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A: matriz com os dados
//Variáveis de saída:
//alfa: vetor das soluções de X^TX*alfa = X^TY
//X: matriz com a primeira coluna de 1's e as demais com os dados
//////////////////////////////////////////////////////////////////////////

function [alfa, X]=Minimos_Quadrados_Manchine_Learning(A)
    [n]=size(A,1);
    X = [ones(n, 1) A(:,1:10)];
    Y = A(:,11);
    [alfa, C, P]=Gaussian_Elimination_4(X'*X,X'*Y);
endfunction

// versão para a questão 3
function [alfa, X]=Minimos_Quadrados_Manchine_Learning_2(A)
    [n]=size(A,1);
    X = [ones(n, 1) A(:,1:5) A(:,8)];
    Y = A(:,11);
    [alfa, C, P]=Gaussian_Elimination_4(X'*X,X'*Y);
endfunction

//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//X: matriz com a primeira coluna de 1's e as demais com os dados
//Y: matriz dos valores a serem previstos
//alfa: vetor obtido pelo métodu dos mínimos quadrados
//Variáveis de saída:
//diag_Prev: diagnóstico previsto (+ ou - 1)
//percent_acerto: porcentagem de acertos em relação à Y
//////////////////////////////////////////////////////////////////////////

function [diag_Prev, percent_acerto]=Classificador(alfa, X, Y)
    [n]=size(Y) //n = 280
    diag_Prev=X*alfa; // aplicando na função h(x)
    diag_Prev=(diag_Prev>=0)*1+(diag_Prev<0)*-1; // substituindo por 1 e -1
    percent_acerto = diag_Prev.*Y // se for igual é 1*1 ou (-1)*(-1), que é 1, e se for diferente é 1*(-1) ou -1*1, que é -1
    percent_acerto = sum(percent_acerto==1)/280;
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
