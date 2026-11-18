#include <cmath>
#include <numbers>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>  

#include "dual.hpp"
#include "exponent.hpp"
#include "flow.hpp"
#include "reference.hpp"

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

TEST_CASE("cos on Dual") {
	const Dual c = cos(Dual(0.5, 1.0));
	CHECK_THAT(c.value, WithinAbs(std::cos(0.5), 1e-15));
	CHECK_THAT(c.derivative, WithinAbs(-std::sin(0.5), 1e-15));
}

TEST_CASE("Onedim maps") {
	CHECK_THAT(map_exponent(Logistic{ 2.5 }, 1000000), WithinAbs(std::log(0.5), 1e-5));
	CHECK_THAT(map_exponent(Quadratic{ 2.0 }, 1000000), WithinAbs(std::log(2.0), 1e-5));
}

TEST_CASE("Generic exponent") {
	CHECK_THAT(map_exponent(Alternating{ 2.5, 3.5 }, 1000000),
		WithinAbs(alternating_exponent(2.5, 3.5, 1000000), 1e-6));
}

TEST_CASE("Henon map") {
	CHECK_THAT(plane_exponent(Henon{ 1.4, 0.3 }, 1000000), WithinAbs(0.419222, 0.001));
	CHECK_THAT(plane_exponent(Henon{ 2.0, 0.0 }, 1000000), WithinAbs(std::log(2.0), 1e-5));
	CHECK(std::isnan(plane_exponent(Henon{ 2.0, 0.3 }, 10000000)));
}

TEST_CASE("Henon spectum") {
	const Spectrum both = plane_spectrum(Henon{ 1.4, 0.3 }, 1000000);
	CHECK_THAT(both.first + both.second, WithinAbs(std::log(0.3), 1e-9));
	CHECK_THAT(1.0 - both.first / both.second, WithinAbs(1.258269, 0.001));
}

TEST_CASE("Runge-Kutta") {
	const Lorenz lorenz{ 10.0, 28.0, 8.0 / 3.0 };
	const Vec3<double> start{ 1.0, 1.0, 1.0 };
	const double coarse = flow(lorenz, start, 0.5, 0.0025).x;
	const double medium = flow(lorenz, start, 0.5, 0.00125).x;
	const double fine = flow(lorenz, start, 0.5, 0.000625).x;
	CHECK_THAT((coarse - medium) / (medium - fine), WithinAbs(16.0, 1.0));
}

TEST_CASE("Lorenz exponents") {
	const Lorenz lorenz{ 10.0, 28.0, 8.0 / 3.0 };
	// https://arxiv.org/pdf/1706.05882 0.9056
	CHECK_THAT(flow_exponent(lorenz, 0.01, 1000000), WithinAbs(0.9056, 0.002));
	const Spectrum3 all = flow_spectrum(lorenz, 0.01, 1000000);
	CHECK_THAT(all.first + all.second + all.third, WithinAbs(-(10.0 + 1.0 + 8.0 / 3.0), 0.001));
	CHECK_THAT(all.second, WithinAbs(0.0, 0.001));
	CHECK_THAT(2.0 - (all.first + all.second) / all.third, WithinAbs(2.06, 0.01));
}

TEST_CASE("Duffing") {
	const Duffing duffing{ 0.3, -1.0, 1.0, 0.5, 1.2 };
	const double period = 2.0 * std::numbers::pi / duffing.omega;
	const Spectrum3 forced = flow_spectrum(duffing, period / 100.0, 1000000);
	CHECK_THAT(forced.first + forced.second, WithinAbs(-0.3, 1e-5));
	CHECK_THAT(forced.third, WithinAbs(0.0, 1e-12));
}
