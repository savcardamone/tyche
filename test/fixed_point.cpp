/**
 * @file fixed_point.cpp
 * @author Salvatore Cardamone
 * @brief Unit test for FixedPoint class.
 */

#define BOOST_TEST_MODULE FixedPointTest
#include <boost/test/included/unit_test.hpp>
#include "utilities/fixed_point.hpp"

BOOST_AUTO_TEST_CASE( constructors ) {

  // Q4.4 representation, unsigned
  tycheplusplus::FixedPoint<unsigned char,4> uchar(15.0);
  BOOST_CHECK_EQUAL(uchar.AsDouble(), 15.0);
  // Q3.4 representation, signed
  tycheplusplus::FixedPoint<signed char,4> char_pos(7.0);
  tycheplusplus::FixedPoint<signed char,4> char_neg(-7.0);
  BOOST_CHECK_EQUAL(char_pos.AsDouble(), 7.0);
  BOOST_CHECK_EQUAL(char_neg.AsDouble(), -7.0);

  
}

BOOST_AUTO_TEST_CASE( addition ) {

  tycheplusplus::FixedPoint<unsigned char,4> a(5.0), b(6.0);
  auto c = a + b;
  BOOST_CHECK_EQUAL(c.AsDouble(), 11.0);

}
