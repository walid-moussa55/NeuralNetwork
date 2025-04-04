#include <iostream>
#include <stdexcept>
#include <vector>
#include <initializer_list>
#include <string>
#include <sstream>


template <typename T>
class Matrix{
private:
    std::vector<T> m;
    int m_h; // number of rows
    int m_l; // number of columns

public:
    Matrix(const int h, const int l)
    :   m_h(h), m_l(l) {
        m.resize(m_h*m_l);
    }

    Matrix(const Matrix& m2)
    :   m_h(m2.get_h()), m_l(m2.get_l()), m(m2.m) {}

    Matrix(const std::initializer_list<T> arr)
    :   m_h(arr.size()), m_l(1) {
        m.resize(m_h*m_l);
        int i = 0;
        for (const auto& element : arr){
            m[i++] = element;
        }
    }

    Matrix(const std::initializer_list<std::initializer_list<T>> arr)
    :   m_h(arr.size()), m_l(arr.begin()->size()) {
        m.resize(m_h*m_l);
        int i = 0;
        for (const auto& row : arr){
            if (row.size() != m_l){
                throw std::invalid_argument("Matrix dimensions must agree");
            }
            for (const auto& element : row){
                m[i++] = element;
            }
        }
    }

    ~Matrix(){
        m.clear();
    }

    T& operator()(int i, int j){
        return m[i*m_l + j];
    }

    T operator()(int i, int j) const{
        return m[i*m_l + j];
    }

    int get_h() const{
        return m_h;
    }

    int get_l() const{
        return m_l;
    }

    std::string shape() const{
        std::stringstream ss;
        ss << "(" << m_h << "," << m_l << ")";
        return ss.str();
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& m){
        os << "[";
        for(int i = 0; i < m.get_h(); i++){
            if (i == 0)
                os << "[ ";
            else
                os << " [ ";
            for(int j = 0; j < m.get_l(); j++){
                os << m(i, j) << " ";
            }
            os << "]";
            if (i < m.get_h() - 1)
                os << std::endl;
        }
        os << "]";
        return os;
    }

    Matrix<T> operator+(const Matrix<T>& m2) const{
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") + (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) + m2(i, j);
            }
        }
        return res;
    }
    Matrix<T> operator-(const Matrix<T>& m2) const{
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") - (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) - m2(i, j);
            }
        }
        return res;
    }
    Matrix<T> operator*(const Matrix<T>& m2) const{
        if(m_l != m2.get_h()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") * (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m_h, m2.get_l());
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m2.get_l(); j++){
                res(i, j) = 0;
                for(int k = 0; k < m_l; k++){
                    res(i, j) += (*this)(i, k) * m2(k, j);
                }
            }
        }
        return res;
    }
    Matrix<T> operator/(const Matrix<T>& m2) const{
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") / (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                if(m2(i, j) == 0){
                    throw std::invalid_argument("Division by zero");
                }
                res(i, j) = (*this)(i, j) / m2(i, j);
            }
        }
        return res;
    }

    Matrix<T> operator+(T a) const{
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) + a;
            }
        }
        return res;
    }
    Matrix<T> operator-(T a) const{
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) - a;
            }
        }
        return res;
    }
    Matrix<T> operator*(T a) const{
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) * a;
            }
        }
        return res;
    }
    Matrix<T> operator/(T a) const{
        if(a == 0){
            throw std::invalid_argument("Division by zero");
        }
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = (*this)(i, j) / a;
            }
        }
        return res;
    }

    Matrix<T>& operator+=(const Matrix<T>& m2){
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") += (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) += m2(i, j);
            }
        }
        return *this;
    }

    Matrix<T>& operator-=(const Matrix<T>& m2){
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") -= (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) -= m2(i, j);
            }
        }
        return *this;
    }

    Matrix<T>& operator*=(const Matrix<T>& m2){
        if(m_l != m2.get_h()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m_h) + "," + std::to_string(m_l) + ") *= (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m_h, m2.get_l());
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m2.get_l(); j++){
                res(i, j) = 0;
                for(int k = 0; k < m_l; k++){
                    res(i, j) += (*this)(i, k) * m2(k, j);
                }
            }
        }
        return *this;
    }
    
    Matrix<T>& operator+=(T a){
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) += a;
            }
        }
        return *this;
    }
    Matrix<T>& operator-=(T a){
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) -= a;
            }
        }
        return *this;
    }
    Matrix<T>& operator*=(T a){
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) *= a;
            }
        }
        return *this;
    }
    Matrix<T>& operator/=(T a){
        if(a == 0){
            throw std::invalid_argument("Division by zero");
        }
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                (*this)(i, j) /= a;
            }
        }
        return *this;
    }

    Matrix<T> operator-() const{
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = -(*this)(i, j);
            }
        }
        return res;
    }

    Matrix<T> operator+() const{
        Matrix<T> res(m_h, m_l);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(i, j) = +(*this)(i, j);
            }
        }
        return res;
    }

    bool operator==(const Matrix<T>& m2) const{
        if(m_h != m2.get_h() || m_l != m2.get_l()){
            return false;
        }
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                if((*this)(i, j) != m2(i, j)){
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const Matrix<T>& m2) const{
        return !(*this == m2);
    }

    Matrix<T>& operator=(const Matrix<T>& m2){
        m_h = m2.get_h();
        m_l = m2.get_l();
        m = m2.m;
        return *this;
    }

    Matrix<T> operator[](int i) const{
        if(i < 0 || i >= m_h){
            throw std::out_of_range("Index out of range");
        }
        Matrix<T> res(1, m_l);
        for(int j = 0; j < m_l; j++){
            res(0, j) = (*this)(i, j);
        }
        return res;
    }

    // transpose matrix
    Matrix<T> transpose() const{
        Matrix<T> res(m_l, m_h);
        for(int i = 0; i < m_h; i++){
            for(int j = 0; j < m_l; j++){
                res(j, i) = (*this)(i, j);
            }
        }
        return res;
    }

    // random matrix
    static Matrix<T> random(const int h, const int l, T min, T max){
        Matrix<T> res(h, l);
        for(int i = 0; i < h; i++){
            for(int j = 0; j < l; j++){
                // Use static_cast to double for floating-point division
                double random_value = static_cast<double>(rand()) / RAND_MAX;
                res(i, j) = min + static_cast<T>(random_value * (max - min));
            }
        }
        return res;
    }

    // identity matrix
    static Matrix<T> identity(const int n){
        Matrix<T> res(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                res(i, j) = (i == j) ? 1 : 0;
            }
        }
        return res;
    }

    // zero matrix
    static Matrix<T> zeros(const int h, const int l){
        Matrix<T> res(h, l);
        for(int i = 0; i < h; i++){
            for(int j = 0; j < l; j++){
                res(i, j) = 0;
            }
        }
        return res;
    }

    // one matrix
    static Matrix<T> ones(const int h, const int l){
        Matrix<T> res(h, l);
        for(int i = 0; i < h; i++){
            for(int j = 0; j < l; j++){
                res(i, j) = 1;
            }
        }
        return res;
    }

    static Matrix<T> multiply(const Matrix<T>& m1, const Matrix<T>& m2){
        if(m1.get_h() != m2.get_h() || m1.get_l() != m2.get_l()){
            throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m1.get_h()) + "," + std::to_string(m1.get_l()) + ") multiply (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
        }
        Matrix<T> res(m1.get_h(), m1.get_l());
        for(int i = 0; i < m1.get_h(); i++){
            for(int j = 0; j < m1.get_l(); j++){
                res(i, j) = m1(i, j) * m2(i, j);
            }
        }
        return res;
    }

    static Matrix<T> determinant(const Matrix<T>& m){
        if(m.get_h() != m.get_l()){
            throw std::invalid_argument("Matrix must be square");
        }
        if(m.get_h() == 1){
            return m;
        }
        if(m.get_h() == 2){
            return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
        }
        T det = 0;
        for(int i = 0; i < m.get_h(); i++){
            Matrix<T> submatrix(m.get_h() - 1, m.get_l() - 1);
            for(int j = 1; j < m.get_h(); j++){
                for(int k = 0; k < m.get_l(); k++){
                    if(k < i){
                        submatrix(j - 1, k) = m(j, k);
                    }
                    else if(k > i){
                        submatrix(j - 1, k - 1) = m(j, k);
                    }
                }
            }
            det += m(0, i) * determinant(submatrix) * ((i % 2 == 0) ? 1 : -1);
        }
        return det;
    }

    static Matrix<T> inverse(const Matrix<T>& m){
        if(m.get_h() != m.get_l()){
            throw std::invalid_argument("Matrix must be square");
        }
        T det = determinant(m);
        if(det == 0){
            throw std::invalid_argument("Matrix is singular");
        }
        Matrix<T> res(m.get_h(), m.get_l());
        for(int i = 0; i < m.get_h(); i++){
            for(int j = 0; j < m.get_l(); j++){
                Matrix<T> submatrix(m.get_h() - 1, m.get_l() - 1);
                for(int k = 0; k < m.get_h(); k++){
                    for(int l = 0; l < m.get_l(); l++){
                        if(k < i && l < j){
                            submatrix(k, l) = m(k, l);
                        }
                        else if(k < i && l > j){
                            submatrix(k, l - 1) = m(k, l);
                        }
                        else if(k > i && l < j){
                            submatrix(k - 1, l) = m(k, l);
                        }
                        else if(k > i && l > j){
                            submatrix(k - 1, l - 1) = m(k, l);
                        }
                    }
                }
                res(j, i) = determinant(submatrix) * (((i + j) % 2 == 0) ? 1 : -1) / det;
            }
        }
        return res;
    }

    static Matrix<T> mean(const Matrix<T>& m){
        Matrix<T> res(1, m.get_l());
        for(int i = 0; i < m.get_l(); i++){
            T sum = 0;
            for(int j = 0; j < m.get_h(); j++){
                sum += m(j, i);
            }
            res(0, i) = sum / m.get_h();
        }
        return res;
    }

    static Matrix<T> sum(const Matrix<T>& m){
        Matrix<T> res(1, m.get_l());
        for(int i = 0; i < m.get_l(); i++){
            T sum = 0;
            for(int j = 0; j < m.get_h(); j++){
                sum += m(j, i);
            }
            res(0, i) = sum;
        }
        return res;
    }

    static Matrix<T> concatenate(const Matrix<T>& m1, const Matrix<T>& m2, int axis){
        if(axis == 0){
            if(m1.get_l() != m2.get_l()){
                throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m1.get_h()) + "," + std::to_string(m1.get_l()) + ") concatenate (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
            }
            Matrix<T> res(m1.get_h() + m2.get_h(), m1.get_l());
            for(int i = 0; i < m1.get_h(); i++){
                for(int j = 0; j < m1.get_l(); j++){
                    res(i, j) = m1(i, j);
                }
            }
            for(int i = 0; i < m2.get_h(); i++){
                for(int j = 0; j < m2.get_l(); j++){
                    res(i + m1.get_h(), j) = m2(i, j);
                }
            }
            return res;
        }
        else if(axis == 1){
            if(m1.get_h() != m2.get_h()){
                throw std::invalid_argument("Matrix dimensions must agree, (" + std::to_string(m1.get_h()) + "," + std::to_string(m1.get_l()) + ") concatenate (" + std::to_string(m2.get_h()) + "," + std::to_string(m2.get_l()) + ")");
            }
            Matrix<T> res(m1.get_h(), m1.get_l() + m2.get_l());
            for(int i = 0; i < m1.get_h(); i++){
                for(int j = 0; j < m1.get_l(); j++){
                    res(i, j) = m1(i, j);
                }
            }
            for(int i = 0; i < m2.get_h(); i++){
                for(int j = 0; j < m2.get_l(); j++){
                    res(i, j + m1.get_l()) = m2(i, j);
                }
            }
            return res;
        }
        else{
            throw std::invalid_argument("Axis must be 0 or 1");
        }
    }
};

template <typename T>
Matrix<T> operator+(T a, const Matrix<T>& m){
    return m + a;
}

template <typename T>
Matrix<T> operator-(T a, const Matrix<T>& m){
    return -m + a;
}

template <typename T>
Matrix<T> operator*(T a, const Matrix<T>& m){
    return m * a;
}
