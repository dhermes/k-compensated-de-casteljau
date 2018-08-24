#include <cmath>
#include <tuple>
#include <vector>

namespace eft {

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
}

namespace de_casteljau {

double basic(double s, std::vector<double> coeffs)
{
    size_t degree, j, k;
    double r;
    std::vector<double> pk(coeffs);

    r = 1.0 - s;
    degree = coeffs.size() - 1;
    for (k = 0; k < degree; ++k) {
        for (j = 0; j < degree - k; ++j) {
            pk[j] = r * pk[j] + s * pk[j + 1];
        }
    }

    return pk[0];
}
}
