/**
 * @file fixed_point.hpp
 * @author Salvatore Cardamone
 * @brief Largely a replica of the Fixed Point Math Library developed by
 *        Peter Schregle. Pared down superfluous functionality and added some
 *        bits and pieces to make the code a little more useful for the types
 *        of calculation we do in tyche++.
 */
#ifndef __TYCHEPLUSPLUS_FIXED_POINT_HPP
#define __TYCHEPLUSPLUS_FIXED_POINT_HPP

#include <boost/operators.hpp>
#include "utilities/type_promotion.hpp"
#include <boost/type_index.hpp>

namespace tycheplusplus {
  
/**
 * @class FixedPoint
 * @brief Generic fixed point functionality, allowing us to use this class as
 *        we would any other data type, like float or double. All integer
 *        arithmetic, so if there's no dedicated FPU, using FixedPoint as the
 *        real numerical type is probably advantageous. Furthermore, any
 *        high-level synthesis tools will hopefully pick up on the use of
 *        integer arithmetic and generate efficient implementations.
 *
 *        boost/operators.hpp provides a fairly remarkable set of
 *        functionalities, whereby operators can be automatically generated from
 *        a smaller set of operators. Deriving from the appropriate boost
 *        classes within boost/operators.hpp then gives us access to the full
 *        complement of operators.
 *
 *        For instance, deriving from boost::ordered_field_operators allows us
 *        to explicitly define the += operator, and obtain the + operator for
 *        free, without the need for additional boilerplate.
 *
 * @tparam B Data type used to store the fixed point number. It the type is
 *           signed, then the fixed point representation will be signed too.
 * @tparam I Number of integer bits.
 * @tparam F Number of fractional bits. Automatically determined from number of
 *           integer bits and number of available storage bits.
 */
template<typename B,
         unsigned char I,
         unsigned char F = std::numeric_limits<B>::digits - I>
class FixedPoint
    : boost::ordered_field_operators<FixedPoint<B,I,F>,
      boost::unit_steppable<FixedPoint<B,I,F>,
      boost::shiftable<FixedPoint<B,I,F> > > >
{

public:
  /**
   * @brief Class constructor.
   * @param value Single precision value to initialise with.
   */
  FixedPoint(float value)
      : value_(value * two_power_f_ + (value >= 0 ? 0.5 : -0.5)) {}

  /**
   * @brief Class constructor.
   * @param value Double precision value to initialise with.
   */
  FixedPoint(double value)
      : value_(value * two_power_f_ + (value >= 0 ? 0.5 : -0.5)) {
  }

  /**
   * @brief Class constructor.
   * @param value Single precision value to initialise with.
   */
  fixed& operator +=(fixed const& rhs) {
    value_ += rhs.value_;
    return *this;
  }

  /**
   * @brief Subtraction assignment operator.
   * @param rhs Value to subtract from lhs.
   * @retval lhs - rhs.
   */
  fixed& operator -=(fixed const& rhs) {
    value_ -= rhs.value_;
    return *this;
  }

  /**
   * @brief Multiplication assignment operator.
   * @param rhs Value to multiply lhs by.
   * @retval lhs * rhs.
   */
  fixed& operator *=(fixed const& rhs) {
    value_ = (static_cast<typename TypePromotion<B>::type>
	      (value_) * rhs.value_) >> number_fractional_bits_;
    return *this;
  }

  /**
   * @brief Division assignment operator.
   * @param rhs Value to divide lhs by.
   * @retval lhs / rhs.
   */
  fixed& operator /=(fixed const& rhs) {
    value_ = (static_cast<typename TypePromotion<B>::type>
	      (value_) << number_fractional_bits_) / rhs.value_;
    return *this;
  }
  
  /**
   * @brief Convert the internal value of the fixed point object to a float.
   * @retval Fixed point number cast to float.
   */
  float AsFloat() const {
    return (float)value_ / two_power_f_ ;
  }

  /**
   * @brief Convert the internal value of the fixed point object to a double.
   * @retval Fixed point number cast to double.
   */
  double AsDouble() const {
    return (double)value_ / two_power_f_ ;
  }

  /**
   * @brief Print some information about the object.
   */
  void Print(std::ostream& stream) const {
    stream << " *** FixedPoint object" << std::endl
	   << "     "
	   << (int)number_fractional_bits_ << " fractional bits and "
	   << (int)number_integer_bits_ << " integer bits." << std::endl
	   << "     Storage type: "
	   << boost::typeindex::type_id<B>().pretty_name() << std::endl
	   << "     Has Sign Bit: "
	   << std::numeric_limits<B>::is_signed << std::endl
	   << "     Stored Value: " << value_ << std::endl
	   << "     Floating Point: " << AsDouble() << std::endl;
  }
  
private:
  // Alias for the sake of simplifying functions
  using fixed = FixedPoint<B,I,F>;
  static constexpr unsigned char number_integer_bits_ = (unsigned char)I;
  static constexpr unsigned char number_fractional_bits_ = (unsigned char)F;
  static constexpr B two_power_f_ = (1ULL << F);
  B value_;

} ;

}

#endif /* #ifndef __TYCHEPLUSPLUS_FIXED_POINT_HPP */
