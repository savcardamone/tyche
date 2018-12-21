/**
 * @file type_promotion.hpp
 * @author Salvatore Cardamone
 * @brief Some template metaprogramming magic that lets us promote data types
 *        to the "next one up" at compile time.
 */
#ifndef __TYCHEPLUSPLUS_TYPE_PROMOTION_HPP
#define __TYCHEPLUSPLUS_TYPE_PROMOTION_HPP

namespace tycheplusplus {

/**
 * @class TypePromotion
 * @brief The functionality of this class is largely equivalent to the boost
 *        type_promotion implementation, but allows us to promote 32-bit types
 *        as well.
 *
 *        Given a particular datatype upon which this class is templated, the
 *        type field contains the datatype which is the next power of two up.
 *        So templating on int will result in the type field containing a 
 *        long long, and so on.
 *
 *        Signed types are promoted to other signed types, and unsigned are
 *        promoted to unsigned. We can promote anything up to and including
 *        64-bit data types, although the 128-bit promotions utilise the
 *        gcc 128-bit datatypes. Apparently clang can also utilise this 
 *        datatype, but I'll need to test.
 *
 *        If no type promotion is supported (floating point, for instance), then
 *        void is returned. Perhaps in future revisions I should learn how to 
 *        generate a compilation error.
 * @ref https://stackoverflow.com/questions/16168927/
 */
template<typename T>
struct TypePromotion {

  typedef
  typename std::conditional<std::is_same<unsigned char,T>::value,unsigned short,
  typename std::conditional<std::is_same<signed char,T>::value,signed short,
  typename std::conditional<std::is_same<unsigned short,T>::value,unsigned long,
  typename std::conditional<std::is_same<signed short,T>::value,signed long,
  typename std::conditional<std::is_same<unsigned int,T>::value,unsigned long int,
  typename std::conditional<std::is_same<int,T>::value,long int,
  typename std::conditional<std::is_same<unsigned long int,T>::value,__uint128_t,
  typename std::conditional<std::is_same<long int,T>::value,__int128_t,
  void>::type>::type>::type>::type>::type>::type>::type>::type type;
			    
};

}

#endif /* #ifndef __TYCHEPLUSPLUS_TYPE_PROMOTION_HPP */
