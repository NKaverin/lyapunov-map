#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>  

#include "dual.hpp"
#include "exponent.hpp"

using Catch::Matchers::WithinAbs;

TEST_CASE("Dual + derivative") {
	const Dual x(0.3, 1.0);
	const Dual y = 4.0 * x * (1.0 - x);
	CHECK_THAT(y.value, WithinAbs(0.84, 1e-12));
	CHECK_THAT(y.derivative, WithinAbs(1.6, 1e-12));
}

TEST_CASE("Logistic map") {
	CHECK_THAT(map_exponent(Logistic{ 4.0 }, 1000000), WithinAbs(std::log(2.0), 1e-5));
}