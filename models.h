#include <math.h>

// Linear Regression Model
// y = W0 + W1*x1 + W2*x2 + ... + Wn*xn
class LinearRegression{
    private:
        Matrix<double> W;
        double learning_rate;
        int eterations;
    public:
        LinearRegression(int n_parameters, double learning_rate, int eterations)
        :   W(Matrix<double>::random(n_parameters+1, 1, 0, 10)),
            learning_rate(learning_rate),
            eterations(eterations) {}
    
        void fit(const Matrix<double>& X, const Matrix<double>& Y){
            Matrix<double> X_aug = Matrix<double>::concatenate(X, Matrix<double>::ones(X.get_h(), 1), 1);
            for(int i = 0; i < eterations; i++){
                Matrix<double> error = X_aug * W - Y;
                Matrix<double> gradient = X_aug.transpose() * error / X_aug.get_h();
                Matrix<double> loss = MSE(Y, X_aug * W);
                std::cout <<"iteration: "<< i+1 << " Loss: " << loss << std::endl;
                W = W - gradient * learning_rate;
            }
        }
    
        Matrix<double> predict(const Matrix<double>& X){
            Matrix<double> X_aug = Matrix<double>::concatenate(X, Matrix<double>::ones(X.get_h(), 1), 1);
            return X_aug * W;
        }
    
        Matrix<double> get_weights(){
            return W;
        }
    private:
        Matrix<double> MSE(const Matrix<double>& Y_true, const Matrix<double>& Y_pred){
            return Matrix<double>::mean((Y_true - Y_pred).transpose() * (Y_true - Y_pred));
        }
};

// Polynomial Regression Model
// y = W0 + W1*x1 + W2*x2 + ... + Wn*xn
class PolynomialRegression{
    private:
        Matrix<double> W;
        double learning_rate;
        int eterations;
        int order;
    public:
        PolynomialRegression(int order, double learning_rate, int eterations)
        : W(Matrix<double>::random(order+1, 1, 0, 10)),
          learning_rate(learning_rate),
          eterations(eterations), order(order) {}

        void fit(const Matrix<double>& X, const Matrix<double>& Y){
            Matrix<double> X_poly = transform(X);
            for(int i = 0; i < eterations; i++){
                Matrix<double> error = X_poly * W - Y;
                Matrix<double> gradient = X_poly.transpose() * error / X_poly.get_h();
                Matrix<double> loss = MSE(Y, X_poly * W);
                std::cout <<"iteration: "<< i+1 << " Loss: " << loss << std::endl;
                W = W - gradient * learning_rate;
            }
        }

        Matrix<double> predict(const Matrix<double>& X){
            Matrix<double> X_poly = transform(X);
            return X_poly * W;
        }

        Matrix<double> get_weights(){
            return W;
        }

    private:
        Matrix<double> transform(const Matrix<double>& X){
            Matrix<double> X_poly(X.get_h(), order+1);
            for(int i = 0; i < X.get_h(); i++){
                for(int j = 0; j < order+1; j++){
                    X_poly(i, j) = std::pow(X(i, 0), j);
                }
            }
            return X_poly;
        }

        Matrix<double> MSE(const Matrix<double>& Y_true, const Matrix<double>& Y_pred){
            return Matrix<double>::mean((Y_true - Y_pred).transpose() * (Y_true - Y_pred));
        }
};

// Artificial Neural Network
class ANN{
    private:
        Matrix<double> w_h;
        Matrix<double> w_o;

    public:
        ANN(int input_size, int hidden_size, int output_size)
        : w_h(Matrix<double>::random(input_size+1, hidden_size, -1, 1)), w_o(Matrix<double>::random(hidden_size+1, output_size, -1, 1))
        {}

        void fit(const Matrix<double>& X_train, const Matrix<double>& y_train, int eterations, double learning_rate){

            Matrix<double> X_aug = Matrix<double>::concatenate(X_train, Matrix<double>::ones(X_train.get_h(), 1), 1);
            // std::cout << "Training..." << std::endl;
            for(int i=0;i<eterations;i++){
                // std::cout << "forwardA" << std::endl;
                Matrix<double> A = forward_A(X_train);
                // std::cout << "forwardY" << std::endl;
                Matrix<double> Y_pred = forward_Y(A);
                // std::cout << "backwardY" << std::endl;
                Matrix<double> A_aug = Matrix<double>::concatenate(A, Matrix<double>::ones(A.get_h(), 1), 1);
                Matrix<double> dw_o = backward_Y(A_aug, y_train, Y_pred);
                // std::cout << "backwardA" << std::endl;
                Matrix<double> dw_h = backward_A(X_aug, y_train, Y_pred);
                // std::cout << "update w_o" << std::endl;
                w_o -= dw_o * learning_rate;
                // std::cout << "update w_h" << std::endl;
                w_h -= dw_h * learning_rate;
                // std::cout << "loss" << std::endl;
                Matrix<double> error = loss(y_train, Y_pred);
                std::cout << "iteration: " << i+1 << " Loss: " << error << std::endl;
            }
        }

        Matrix<double> predict(const Matrix<double>& X){
            // std::cout << "Predicting..." << std::endl;
            Matrix<double> A = forward_A(X);
            Matrix<double> Y_pred = forward_Y(A);
            return Y_pred;
        }

    private:
        Matrix<double> _ReLU(const Matrix<double>& X){
            Matrix<double> res(X.get_h(), X.get_l());
            for(int i=0; i<X.get_h();i++){
                for(int j=0;j<X.get_l();j++){
                    if (X(i,j) > 0){
                        res(i, j) = X(i,j);
                    }else{
                        res(i,j) = 0;
                    }
                }
            }
            return res;
        }

        Matrix<double> forward_A(const Matrix<double>& X){
            // std::cout << "X" << X.shape() << std::endl;
            Matrix<double> X_aug = Matrix<double>::concatenate(X, Matrix<double>::ones(X.get_h(), 1),1);
            // std::cout << "X_aug" << X_aug.shape() << std::endl;
            // std::cout << "w_h" << w_h.shape() << std::endl;
            Matrix<double> H = X_aug * w_h;
            // std::cout << "H" << H.shape() << std::endl;
            Matrix<double> A = _ReLU(H);
            // std::cout << "A" << A.shape() << std::endl;
            return A;
        }

        Matrix<double> forward_Y(const Matrix<double>& A){
            // std::cout << "A" << A.shape() << std::endl;
            Matrix<double> A_aug = Matrix<double>::concatenate(A, Matrix<double>::ones(A.get_h(),1),1);
            // std::cout << "A_aug" << A_aug.shape() << std::endl;
            // std::cout << "w_o" << w_o.shape() << std::endl;
            Matrix<double> H = A_aug * w_o;
            // std::cout << "H" << H.shape() << std::endl;
            return H;
        }

        Matrix<double> loss(const Matrix<double>& Y_pred, const Matrix<double>& Y_true){
            // std::cout << "Y_pred" << Y_pred.shape() << std::endl;
            // std::cout << "Y_true" << Y_true.shape() << std::endl;
            return Matrix<double>::mean((Y_pred - Y_true).transpose() * (Y_pred - Y_true));
        }

        Matrix<double> _sign(const Matrix<double>& X){
            Matrix<double> res(X.get_h(), X.get_l());
            for(int i=0;i<X.get_h();i++){
                for(int j=0;j<X.get_l();j++){
                    res(i,j) = (X(i,j) > 0) ? 1 : 0;
                }
            }
            return res;
        }

        Matrix<double> backward_A(const Matrix<double>& X, const Matrix<double>& Y_true, const Matrix<double>& Y_pred){
            // std::cout << "X" << X.shape() << std::endl;
            // std::cout << "Y_true" << Y_true.shape() << std::endl;
            // std::cout << "Y_pred" << Y_pred.shape() << std::endl;
            // std::cout << "w_h" << w_h.shape() << std::endl;
            // std::cout << "w_o" << w_o.shape() << std::endl;
            int m = X.get_h();
            Matrix<double> w_o_sans_biais(w_o.get_h()-1, w_o.get_l());
            for(int i=0;i<w_o.get_h()-1;i++){
                for(int j=0;j<w_o.get_l();j++){
                    w_o_sans_biais(i,j) = w_o(i,j);
                }
            }
            // std::cout << "w_o_sans_biais" << w_o_sans_biais.shape() << std::endl;
            Matrix<double> delta_1 = Matrix<double>::multiply((Y_pred - Y_true) * w_o_sans_biais.transpose(),_sign(X * w_h));
            // std::cout << "delta_1" << delta_1.shape() << std::endl;
            Matrix<double> dw_h = X.transpose() * delta_1 * 2 / m;
            // std::cout << "dw_h" << dw_h.shape() << std::endl;
            return dw_h;
        }

        Matrix<double> backward_Y(const Matrix<double>& A, const Matrix<double>& Y_true, const Matrix<double>& Y_pred){
            // std::cout << "A" << A.shape() << std::endl;
            // std::cout << "Y_true" << Y_true.shape() << std::endl;
            // std::cout << "Y_pred" << Y_pred.shape() << std::endl;
            int m = A.get_h();
            Matrix<double> delta_2 = Y_pred - Y_true;
            // std::cout << "delta_2" << delta_2.shape() << std::endl;
            Matrix<double> dw_o = A.transpose() * delta_2 * 2 / m;
            // std::cout << "dw_o" << dw_o.shape() << std::endl;
            return dw_o;
        }

};






