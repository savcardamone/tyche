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

#include <cstddef>
#include <boost/operators.hpp>
#include <boost/type_index.hpp>
#include "utilities/type_promotion.hpp"

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
  FixedPoint<B,I,F>& operator +=(FixedPoint<B,I,F> const& rhs) {
    value_ += rhs.value_;
    return *this;
  }

  /**
   * @brief Subtraction assignment operator.
   * @param rhs Value to subtract from lhs.
   * @retval lhs - rhs.
   */
  FixedPoint<B,I,F>& operator -=(FixedPoint<B,I,F> const& rhs) {
    value_ -= rhs.value_;
    return *this;
  }

  /**
   * @brief Multiplication assignment operator.
   * @param rhs Value to multiply lhs by.
   * @retval lhs * rhs.
   */
  FixedPoint<B,I,F>& operator *=(FixedPoint<B,I,F> const& rhs) {
    value_ = (static_cast<typename TypePromotion<B>::type>
	      (value_) * rhs.value_) >> number_fractional_bits_;
    return *this;
  }

  /**
   * @brief Division assignment operator.
   * @param rhs Value to divide lhs by.
   * @retval lhs / rhs.
   */
  FixedPoint<B,I,F>& operator /=(FixedPoint<B,I,F> const& rhs) {
    value_ = (static_cast<typename TypePromotion<B>::type>
	      (value_) << number_fractional_bits_) / rhs.value_;
    return *this;
  }
  
  /**
   * @brief Convert the internal value of the fixed point object to a float.
   * @retval Fixed point number cast to float.
   */
  float AsFloat() const {
    return (float)value_ / two_power_f_;
  }

  /**
   * @brief Convert the internal value of the fixed point object to a double.
   * @retval Fixed point number cast to double.
   */
  double AsDouble() const {
    return (double)value_ / two_power_f_;
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
	   << "     Stored Value: " << std::hex << value_ << std::dec << std::endl
	   << "     Floating Point: " << AsDouble() << std::endl;
  }

  /**
   * @brief Compute the exponential of a fixed point number.
   *
   *        This is fairly inefficient, utilising (I+F) integer multiplications.
   *        The exponential is split into its integer and fractional parts:
   *
   *                         exp(i.f) = exp(i) * exp(f)
   *
   *        For the fractional part, we move down the fractional bits of the
   *        argument, lookup the associated value for the exponential of the
   *        fractional bit and multiply-accumulate if the bit is high, otherwise
   *        it makes no contribution. So, for instance, e^{0.625} is equal to:
   *
   *                    1*exp(0.5) * 0*exp(0.25) * 1*exp(0.125)
   *
   *        the values of the exponential for which are already tabulated.
   *        We compute the integer part in a similar fashion using the
   *        integer part lookup table. If the argument is negative, we
   *        divide-accumulate rather than multiply-accumulate.
   *
   *        We should implement a specialised Gaussian function, since the 
   *        integer part is only relevant over a much smaller dynamic range.
   * @param arg The argument of the exponential function.
   * @retval exp(arg).
   */
  friend FixedPoint<B,I,F> exp(FixedPoint<B,I,F> const& arg) {

    FixedPoint<B,I,F> result(1.0);

    // We start from the MSB in the fractional part and work our way down the
    // number of fractional bits
    for (int i_frac = F-1; i_frac >= 0; --i_frac) {
      if (arg.value_ & 1ULL<<i_frac) {
	result.value_ =
	  (static_cast<typename TypePromotion<B>::type>(result.value_) *
	   (exp_frac_lut[F-i_frac-1] >> (32-F))) >> F;
      }
    }

    // Need to find out whether we're dividing or multiplying for the
    // integer part
    bool is_negative =
      std::numeric_limits<B>::is_signed && ((1ULL << (I+F-1)) & arg.value_);

    // If the number is negative, we need to do some two's complement to get the
    // integer part then work our way up from the LSB
    if (is_negative) {
      B integer_part = ~(arg.value_ >> F) + 1;
      for (int i_int = 0; i_int<I; ++i_int) {
	if (integer_part & 1ULL<<i_int) {
	  result.value_ =
	    (static_cast<typename TypePromotion<B>::type>(result.value_) << F) /
	    (exp_int_lut[i_int] >> (32-F));
	}
      }
    // If the number is positive, we start from the MSB in the integer part and
    // work our way up the number of integer bits
    } else {
      for (int i_int = F; i_int<(I+F); ++i_int) {
	if (arg.value_ & 1ULL<<i_int) {
	  result.value_ =
	    (static_cast<typename TypePromotion<B>::type>(result.value_) *
	     (exp_int_lut[i_int-F]) >> (32-F)) >> F;
	}
      }
    }
    
    return result;
    
  }
  
private:
  // Alias for the sake of simplifying functions
  B value_;
  static constexpr unsigned char number_integer_bits_ = (unsigned char)I;
  static constexpr unsigned char number_fractional_bits_ = (unsigned char)F;
  static constexpr B two_power_f_ = (1ULL << F);

  // exp[0.5], exp[0.25], exp[0.125], etc... in Q32.32
  static constexpr unsigned long exp_frac_lut[32] = {
    0x00000001a61298e2, 0x0000000148b5e3c4, 0x000000012216045b,
    0x000000011082b578, 0x0000000108205601, 0x0000000104080ab5,
    0x0000000102020156, 0x000000010100802b, 0x0000000100802005,
    0x0000000100400801, 0x0000000100200200, 0x0000000100100080,
    0x0000000100080020, 0x0000000100040008, 0x0000000100020002,
    0x0000000100010001, 0x0000000100008000, 0x0000000100004000,
    0x0000000100002000, 0x0000000100001000, 0x0000000100000800,
    0x0000000100000400, 0x0000000100000200, 0x0000000100000100,
    0x0000000100000080, 0x0000000100000040, 0x0000000100000020,
    0x0000000100000010, 0x0000000100000008, 0x0000000100000004,
    0x0000000100000002, 0x0000000100000001
  };
  // exp[1], exp[2], exp[4], etc... in Q32.32
  static constexpr unsigned long exp_int_lut[32] = {
    0x00000002b7e15163, 0x0000000763992e35, 0x0000003699205c4e,
    0x00000ba4f53ea386, 0x0087975e85400100, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff
  };

};
  
}

#endif /* #ifndef __TYCHEPLUSPLUS_FIXED_POINT_HPP */


