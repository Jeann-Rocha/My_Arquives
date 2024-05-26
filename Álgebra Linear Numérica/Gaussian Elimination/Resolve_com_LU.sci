//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//X: solução do sistema AX=B (assumimos que tal solução existe).
//P: matriz de permutação utilizada no escalonamento.
//C: Seja PA=LU a decomposição LU de PA.
//Então C(i,j)=L(i,j) para i>j e C(i,j)=U(i,j) para j>=i.
//////////////////////////////////////////////////////////////////////////
function [X]=Resolve_com_LU(C, P, B)
    [m,n]=size(B);
    L=tril(C, -1) + eye(n,n); // matriz triangular inferior
    U=triu(C); // matriz triangular superior
    
    //Permutando (AX = B => LUX = PAX = PB)
    B=P*B;
    
    //Resolvendo LY = B
    Y=zeros(m,n);
    Y(1,:) = B(1,:);
    for i=2:n
        Y(i,:)=B(i,:)-L(i,1:i-1)*Y(1:i-1,:);
    end
    
    //Resolvendo UX = B
    X=zeros(m,n);
    X(n,:)=Y(n,:)/U(n,n);
    for i=n-1:-1:1
        X(i,:)=(Y(i,:)-U(i,i+1:n)*X(i+1:n,:))/U(i,i);
    end
    
endfunction
