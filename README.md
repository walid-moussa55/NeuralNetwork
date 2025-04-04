# Neural Network Project

This project implements a neural network and other machine learning models in C++. Below is a description of the files in this project and instructions on how to work with them.

## File Descriptions

### 1. `main.cpp`
- **Purpose**: The entry point of the application. It demonstrates the usage of the `ANN` class by generating synthetic data, training the model, and making predictions.
- **Key Features**:
  - Generates random data for training and testing.
  - Trains an artificial neural network (ANN) using the `ANN` class.
  - Outputs predictions and compares them with expected results.

### 2. `models.h`
- **Purpose**: Contains implementations of various machine learning models, including:
  - **Linear Regression**: A simple linear model for regression tasks.
  - **Polynomial Regression**: Extends linear regression to handle polynomial features.
  - **Artificial Neural Network (ANN)**: A feedforward neural network with one hidden layer.
- **Key Features**:
  - Each model includes methods for training (`fit`), prediction (`predict`), and loss calculation.
  - The ANN class supports ReLU activation and backpropagation.

### 3. `matrix.h`
- **Purpose**: Implements a generic `Matrix` class to handle mathematical operations required for machine learning models.
- **Key Features**:
  - Supports matrix addition, subtraction, multiplication, division, and scalar operations.
  - Includes utility methods like transpose, random matrix generation, concatenation, and more.
  - Provides static methods for common operations like mean, sum, determinant, and inverse.

### 4. `run.bat`
- **Purpose**: A batch script to compile and run the project on Windows.
- **Key Features**:
  - Compiles `main.cpp` using `g++`.
  - Runs the resulting executable (`app.exe`).

## How to Work with the Project

### Prerequisites
- A C++ compiler (e.g., GCC or MSVC).
- Windows Command Prompt or a similar terminal.

### Steps to Run
1. **Compile and Run**:
   - Double-click `run.bat` or execute it in the terminal:
     ```bat
     run.bat
     ```
   - This will compile `main.cpp` and execute the resulting program.

2. **Modify the Code**:
   - Update `main.cpp` to test different models or datasets.
   - Adjust hyperparameters (e.g., learning rate, iterations) in the `fit` methods of the models in `models.h`.

3. **Extend the Project**:
   - Add new models or activation functions in `models.h`.
   - Implement additional matrix operations in `matrix.h`.

### Notes
- Ensure the `Matrix` class in `matrix.h` is used for all mathematical operations to maintain consistency.
- The project currently uses synthetic data. Replace it with real datasets by modifying `main.cpp`.

Feel free to explore and enhance the project as needed!
