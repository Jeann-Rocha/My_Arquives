//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//xk: solução do sistema encontrada pelo métodu de Jacobi.
//dif_norm: norma da diferença entre as duas últimas aproximações.
//k: número de iterações efetuadas.
//rk_norm: norma do resíduo.
//Então dif_norm = ||xk-x(k-1))|| e rk_norm = ||rk|| = ||b-Axk||.
//////////////////////////////////////////////////////////////////////////
function [xk, dif_norm, k, rk_norm]=GaussSeidel1(A, b, x0, E, M, type_norm)
    L=tril(A); // triangular inferior de A
    U=triu(A,1); // triangular superior de A, exceto a diagonal
    
    k=1; // variável de iteração
    xk=inv(L)*(-U*x0+b); // vetor da primeira solução aproximada
    
    // realizando as iterações até estar abaixo da tolerância ou até o máximo permitido
    while norm(xk-x0,type_norm)>=E && k<M
        x0=xk; // vetor da (k-1)-ésima solução aproximada armazenado em x0
        xk=inv(L)*(-U*x0+b); // vetor da k-ésima solução aproximada
        k=k+1;
    end
    
    dif_norm=norm(xk-x0,type_norm);
    rk_norm=norm(b-A*xk);
endfunction
