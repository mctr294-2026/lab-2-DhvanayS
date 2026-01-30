#include "roots.hpp" //
#include <cmath>

const int MAX_ITER = 1e6;
const double TOL = 1e-6; //tolerance level higher than test values

bool bisection(std::function<double(double)> f, double a, double b,double *root) { //define function and its arguments
    // bool means the function returns a boolean(true or false)
    // std:: function is a variable that can store a function (for easy storage and passing around)
    // <double(double)> means the function takes in a double and returns a double
    // <int(double,double)> means the functions TAKES IN two doubles (brackets is always input) and returns an int < <> means OUTPUT>
    // *root is a pointer to a double where the found root will be stored

    // evauate functiion at endpoints
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) {
        return false; // no root found if f(a) and f(b) have same sign
    }
    for (int i = 0; i < MAX_ITER; ++i) {
        double c = 0.5 * (a + b); // midpoint
        double fc = f(c); // evaluate function at midpoint
        if (std::fabs(fc) < TOL || (b - a) / 2 < TOL) {
            //std::fabs is absolute value function from cmath
            // || means "or"
            // we check (b-a)/2 < tol to ensure inte+-

            *root = c; // store root in pointer aka put the answer in the variable defined in main.cpp
            // pointers are used to allow functions to modify variables outside their own file/scope
            return true;
        }
        if(fa * fc < 0){
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return false; // no root found within max iterations
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root){
    // you need a double in front of root* because it's a pointer to a double //*variable is the correct syntax to reference actual data stored at the pointer address

    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) { // should give negative value
        return false; // no root found if f(a) and f(b) have same sign
    }
    for (int i = 0; i < MAX_ITER; ++i) {
        double c = a - ((f(
            a)*(b-a))/(fb - fa)); // formula given
        double fc = f(c);
        if (std::fabs(fc) < TOL) { //check if function value at c is close enough to zero
            *root = c; //set root to the x-value of the root, c
            // remember * goes in front of pointer variable
            return true;
        }
        if (fa * fc < 0) { // set b to c if root is in [a,c]
            b = c;
            fb = fc;
        } else { // set a to c if root is in [c,b]
            a = c;
            fa = fc;
        } // then this repeats the for loop
    }
    return false; // no root found within max iterations
}

// N = log_2((b-a)/TOL) = number of iterations needed to achieve tolerance-appropiate root

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double x0, double *root){
    //calls two functions poly1 (FUNCTION f) and poly1_deriv (FUNCTION g)
    // double x0 is initial guess which we have set to -1.0 in main.cpp
    // a is set to -200, b is set to 300

    double x_n = x0; // set initial guess
    if (f(a) * f(b) > 0){ //check intervals are on opposite sides of axis and that our intervals arent at a root
        return false;
    }
    for(int i = 0; i < MAX_ITER; i++){
        double fx_n =f(x_n);
        double gx_n = g(x_n);
        if (std::fabs(g(x_n)) < TOL) { // check derivative isnt too small/zero (to avoid division by zero error)
            return false;
        }
        
        double x_new = x_n - (f(x_n)/g(x_n)); // formula for newton-raphson

        if (x_new < a || x_new > b){ // check new guess is within initial interval
            return false;
        }
        double fx = f(x_new);
        if(std::fabs(fx) < TOL){
            *root = x_new;
            return true;
        }

        x_n = x_new; // update guess for next iteration
        // you don't change the bounds as you did for first 2 methods
    }
    return false; // goes to here if max_iterations are complete and returns false for no root found
}

bool secant(std::function <double(double)> f, double a, double b, double x1, double *root){ // a is taken as the second inital guess
    double x_old = x1+TOL; // first initial guess
    double x_current = x1; // second initial guess

    if(f(a) * f(b) > 0){
        return false;
    }
    for(int i = 0; i < MAX_ITER; i++){
        if (std::fabs((f(x_current) - f(x_old))) < TOL){ // avoid division by zero
            return false; //faster method but downside is that this method can fail
        }
        double x_new = x_current - f(x_current)*((x_current - x_old)/(f(x_current) - f(x_old))); // secant formula
        if (std::fabs(f(x_new)) < TOL){
            *root = x_new;
            return true;
        }
        if(x_new < a || x_new > b){ // check new guess is within initial interval
            return false;
        }
        x_old = x_current;
        x_current = x_new;
    }
    return false; // no root found within max iterations
} 