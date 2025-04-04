#include <iostream>
#include "matrix.h"
#include "models.h"

int main(){
    int n_samples = 100;
    int n_parameters = 1;
    
    Matrix<double> X = Matrix<double>::random(n_samples, n_parameters, 0, 10);
    Matrix<double> Y = 2.0 * X - 10.0 + Matrix<double>::random(n_samples, n_parameters, -1, 1) * 0.1;
    
    ANN model(1, 10, 1);
    model.fit(X, Y, 1000, 0.003);

    Matrix<double> X_test = Matrix<double>::random(4, n_parameters, 0, 10);
    Matrix<double> Y_test = 2.0 * X_test - 10.0;
    Matrix<double> Y_pred = model.predict(X_test);
    std::cout << Matrix<double>::concatenate(X_test, Matrix<double>::concatenate(Y_test, Y_pred,1),1) << std::endl;

    
    return 0;
}