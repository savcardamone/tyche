/**
 * @file fixed_point.cpp
 * @author Salvatore Cardamone
 * @brief Unit test for FixedPoint class.
 */
#define BOOST_TEST_MODULE FixedPointTest
#include <random>
#include <boost/test/included/unit_test.hpp>
#include "utilities/fixed_point.hpp"

/**
 * This case tests construction using all datatypes used for storage of the
 * fixed point representation, and simultaneously minimum and maximum values for
 * the associated representation.
 */
BOOST_AUTO_TEST_CASE( constructors ) {

  // Q4.4 representation, unsigned
  tycheplusplus::FixedPoint<unsigned char,4> us_char(15.0);
  BOOST_CHECK_EQUAL(us_char.AsDouble(), 15.0);
  // Q3.4 representation, signed
  tycheplusplus::FixedPoint<signed char,4> char_pos( 7.0);
  tycheplusplus::FixedPoint<signed char,4> char_neg(-7.0);
  BOOST_CHECK_EQUAL(char_pos.AsDouble(),  7.0);
  BOOST_CHECK_EQUAL(char_neg.AsDouble(), -7.0);

  // Q8.8 representation, unsigned
  tycheplusplus::FixedPoint<unsigned short,8> us_short(255.0);
  BOOST_CHECK_EQUAL(us_short.AsDouble(), 255.0);
  // Q7.8 representation, signed
  tycheplusplus::FixedPoint<signed short,8> short_pos( 127.0);
  tycheplusplus::FixedPoint<signed short,8> short_neg(-127.0);
  BOOST_CHECK_EQUAL(short_pos.AsDouble(),  127.0);
  BOOST_CHECK_EQUAL(short_neg.AsDouble(), -127.0);

  // Q16.16 representation, unsigned
  tycheplusplus::FixedPoint<unsigned int,16> us_int(65535.0);
  BOOST_CHECK_EQUAL(us_int.AsDouble(), 65535.0);
  // Q15.16 representation, signed
  tycheplusplus::FixedPoint<signed int,16> int_pos( 32767.0);
  tycheplusplus::FixedPoint<signed int,16> int_neg(-32767.0);
  BOOST_CHECK_EQUAL(int_pos.AsDouble(),  32767.0);
  BOOST_CHECK_EQUAL(int_neg.AsDouble(), -32767.0);

  // Q32.32 representation, unsigned
  tycheplusplus::FixedPoint<unsigned long int,32> us_long(4294967296.0);
  BOOST_CHECK_EQUAL(us_long.AsDouble(), 4294967296.0);
  // Q31.32 representation, signed
  tycheplusplus::FixedPoint<signed long int,32> long_pos( 2147483647.0);
  tycheplusplus::FixedPoint<signed long int,32> long_neg(-2147483647.0);
  BOOST_CHECK_EQUAL(long_pos.AsDouble(),  2147483647.0);
  BOOST_CHECK_EQUAL(long_neg.AsDouble(), -2147483647.0);

}

BOOST_AUTO_TEST_CASE( addition ) {

  tycheplusplus::FixedPoint<unsigned char,4> a(5.0), b(6.0);
  auto c = a + b;
  BOOST_CHECK_EQUAL(c.AsDouble(), 11.0);

}

/**
 * Test the exponential functionality, both signed and unsigned with positive
 * and negative arguments.
 */
BOOST_AUTO_TEST_CASE( exponential ) {

  constexpr int n_samples = 2048;

  // We'll use the Mersenne Twister for our PRNG
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> us_char_dist(0.0, 1.0);
  for (auto i = 0; i < n_samples; ++i) {
    auto val = us_char_dist(gen);
    tycheplusplus::FixedPoint<unsigned char,4,4> arg(val);
    auto exp_result = exp(arg);
    auto diff = fabs(exp_result.AsDouble() - exp(val));
    BOOST_CHECK_SMALL(diff, 1.0);
  }

  std::uniform_real_distribution<double> us_short_dist(0.0, 2.0);
  for (auto i = 0; i < n_samples; ++i) {
    auto val = us_short_dist(gen);
    tycheplusplus::FixedPoint<unsigned short,8,8> arg(val);
    auto exp_result = exp(arg);
    auto diff = fabs(exp_result.AsDouble() - exp(val));
    BOOST_CHECK_SMALL(diff, 0.2);
  }

  std::uniform_real_distribution<double> us_int_dist(0.0, 3.0);
  for (auto i = 0; i < n_samples; ++i) {
    auto val = us_int_dist(gen);
    tycheplusplus::FixedPoint<unsigned int,16,16> arg(val);
    auto exp_result = exp(arg);
    auto diff = fabs(exp_result.AsDouble() - exp(val));
    BOOST_CHECK_SMALL(diff, 0.01);
  }

  std::uniform_real_distribution<double> us_long_dist(0.0, 4.0);
  for (auto i = 0; i < n_samples; ++i) {
    auto val = us_long_dist(gen);
    tycheplusplus::FixedPoint<unsigned long,32,32> arg(val);
    auto exp_result = exp(arg);
    auto diff = fabs(exp_result.AsDouble() - exp(val));
    BOOST_CHECK_SMALL(diff, 0.0000001);
  }


}
