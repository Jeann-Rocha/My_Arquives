//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//A_test: vetor de dados teste
//A_prev: vetor de dados previstos
//Variáveis de saída:
//M: matriz de confusão (na forma abaixo)
//               A_test
//           |  | 1|-1|
//M -> A_prev| 1| a| b|
//           |-1| c| d|
//////////////////////////////////////////////////////////////////////////

function [M]=Confusin_Matrix(A_test, A_prev)
    M=zeros(2,2);
    [n]=size(A_test,1);
    for i=1:n
        if (A_test(i)==1&&A_prev(i)==1)
            M(1,1)=M(1,1)+1;
        elseif (A_test(i)==1&&A_prev(i)==-1)
            M(2,1)=M(2,1)+1;
        elseif (A_test(i)==-1&&A_prev(i)==1)
            M(1,2)=M(1,2)+1;
        else
            M(2,2)=M(2,2)+1;
        end
    end
endfunction

//////////////////////////////////////////////////////////////////////////
//Variáveis de entrada:
//M: matriz de confusão
//////////////////////////////////////////////////////////////////////////

function Metrics(M)
    TP = M(1,1); // True Positives
    TN = M(2,2); // True Negatives
    FP = M(1,2); // False Positives
    FN = M(2,1); // False Negatives
    
    // métricas 
    accuracy=(TP+TN)/sum(M); // acurácia
    precision=TP/(TP+FP); // precisão
    recall=TP/(TP+FN); // recall
    f1_score=2*(precision*recall)/(precision+recall); //f1 score
    probability_alarm_false=FP/(TP+FP) //probabilidade de falso alarme
    probability_alarm_false_omiss=FN/(TN+FN) //probabilidade de falsa omissão de alarme
    
    disp("Accuracy: " + string(accuracy));
    disp("Precision: " + string(precision));
    disp("Recall: " + string(recall));
    disp("F1 Score: " + string(f1_score));
    disp("Probability False Alarm: " + string(probability_alarm_false));
    disp("Probability Alarm False Omiss: " + string(probability_alarm_false_omiss));
endfunction
