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
