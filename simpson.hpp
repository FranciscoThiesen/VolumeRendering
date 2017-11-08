double DoubleSimpson(double a, double b, double (*f) (double x), double* v);

double AdaptiveSimpson(double a, double b, double (*f) (double x), double tol);