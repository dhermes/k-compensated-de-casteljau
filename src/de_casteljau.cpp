#include <cmath>
#include <tuple>

std::tuple<double, double> two_prod(double val1, double val2)
{
    double product, error;
    product = val1 * val2;
    error = std::fma(val1, val2, -product);
    return std::make_tuple(product, error);
}

std::tuple<double, double> two_sum(double val1, double val2)
{
    double sum, almost_val2, error;
    sum = val1 + val2;
    almost_val2 = sum - val1;
    error = (val1 - (sum - almost_val2)) + (val2 - almost_val2);
    return std::make_tuple(sum, error);
}
