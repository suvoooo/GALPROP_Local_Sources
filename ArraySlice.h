#ifndef ARRAYSLICE_H
#define ARRAYSLICE_H
/** \class ArraySlice
 * \brief A class to slice a normal C array
 *
 * Includes some nice helper functions for its size and mathematical
 * operations.
 */

#include <iostream>
#include <stdexcept>
#include <valarray>

template<typename T>
class ArraySlice {
   public:
      /** \brief Constructor taking a pointer and a size
       *
       * The pointer is supposed to be a part of a c array that extends at
       * least to pointer+size.
       */
      ArraySlice ( T *pointer, size_t size ) : fbegin(pointer), fsize(size), fbeginConst(pointer) {}
      ArraySlice ( const T *pointer, size_t size ) :fbegin(NULL), fsize(size), fbeginConst(pointer) {}

      //! Copy constructor, copy the pointer and the size
      ArraySlice (const ArraySlice<T> &old) : fsize(old.fsize), fbegin(old.fbegin), fbeginConst(old.fbeginConst) {}

      //! Return the size of the slice
      size_t size() const { return fsize; }
      
      /** \brief Return a reference to element number i in the Slice
       * 
       * That is element pointer+i in the original
       */ 
      T & operator [] (size_t i) { return *(fbegin+i); }
      const T & operator [] (size_t i) const { return *(fbeginConst+i); }

      //! Sum all the elements
      T sum() {
	 T sum = 0;
	 for (size_t i = 0; i < fsize; ++i) {
	    sum += *(fbegin+i);
	 }
	 return sum;
      }

      //Assignment operators, work by assigning values to each element
      template <typename C>
	 ArraySlice<T> & operator = (C number) {
	    for (size_t i = 0; i < fsize; ++i) {
	       *(fbegin+i) = number;
	    }
	    return (*this);
	 }
      //We need a special assignment operator for same type slices, otherwise
      //the useless default is used
      ArraySlice<T> & operator = (const ArraySlice<T> & slice) {
	 if (fsize != slice.fsize) {
	    std::string errMsg = "Slices need to have same size in assignment";
	    std::cerr<<errMsg<<std::endl;
	    throw(std::invalid_argument(errMsg));
	 }
	 for (size_t i = 0; i < fsize; ++i) {
	    *(fbegin+i) = slice[i];
	 }
	 return (*this);
      }
      template <typename C>
	 ArraySlice<T> & operator = (const ArraySlice<C> & slice) {
	    if (fsize != slice.fsize) {
	       std::string errMsg = "Slices need to have same size in assignment";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i) {
	       *(fbegin+i) = slice[i];
	    }
	    return (*this);
	 }
      template <typename C>
	 ArraySlice<T> & operator = (const std::valarray<C> & valarray) {
	    if (fsize != valarray.size()) {
	       std::string errMsg = "Valarray and slice need to have same size in assignment";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) = valarray[i];
	    }
	    return (*this);
	 }


      //Mathematical operators should be self-explanatory
      //Includes operation with valarrays and regular numbers
      
      //Addition with a single number
      template <typename C>
	 ArraySlice<T> & operator += (C number) {
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) += number;
	    }
	    return (*this);
	 }
      //Subtraction with a single number
      template <typename C>
	 ArraySlice<T> & operator -= (C number) {
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) -= number;
	    }
	    return (*this);
	 }
      //Multiplication with a single number
      template <typename C>
	 ArraySlice<T> & operator *= (C number) {
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) *= number;
	    }
	    return (*this);
	 }
      //Division with a single number
      template <typename C>
	 ArraySlice<T> & operator /= (C number) {
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) /= number;
	    }
	    return (*this);
	 }

      //Addition with another slice
      template <typename C>
	 ArraySlice<T> & operator += (const ArraySlice<C> & slice) {
	    if (fsize != slice.fsize) {
	       std::string errMsg = "Slices need to have same size in addition";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) += slice[i];
	    }
	    return (*this);
	 }
      //Subtraction with another slice
      template <typename C>
	 ArraySlice<T> & operator -= (const ArraySlice<C> & slice) {
	    if (fsize != slice.fsize) {
	       std::string errMsg = "Slices need to have same size in subtraction";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) -= slice[i];
	    }
	    return (*this);
	 }
      //Multiplication with another slice
      template <typename C>
	 ArraySlice<T> & operator *= (const ArraySlice<C> & slice) {
	    if (fsize != slice.fsize) {
	       std::string errMsg = "Slices need to have same size in multiplication";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) *= slice[i];
	    }
	    return (*this);
	 }
      //Division with another slice
      template <typename C>
	 ArraySlice<T> & operator /= (const ArraySlice<C> & slice) {
	    if (fsize != slice.fsize) {
	       std::string errMsg = "Slices need to have same size in division";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) /= slice[i];
	    }
	    return (*this);
	 }

      //Addition with another valarray
      template <typename C>
	 ArraySlice<T> & operator += (const std::valarray<C> & valarray) {
	    if (fsize != valarray.size()) {
	       std::string errMsg = "Valarray and slice need to have same size in addition";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) += valarray[i];
	    }
	    return (*this);
	 }
      //Subtraction with another valarray
      template <typename C>
	 ArraySlice<T> & operator -= (const std::valarray<C> & valarray) {
	    if (fsize != valarray.size()) {
	       std::string errMsg = "Valarray and slice need to have same size in subtraction";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) -= valarray[i];
	    }
	    return (*this);
	 }
      //Multiplication with another valarray
      template <typename C>
	 ArraySlice<T> & operator *= (const std::valarray<C> & valarray) {
	    if (fsize != valarray.size()) {
	       std::string errMsg = "Valarry and slice need to have same size in multiplication";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) *= valarray[i];
	    }
	    return (*this);
	 }
      //Division with another valarray
      template <typename C>
	 ArraySlice<T> & operator /= (const std::valarray<C> & valarray) {
	    if (fsize != valarray.size()) {
	       std::string errMsg = "Valarray and slice need to have same size in division";
	       std::cerr<<errMsg<<std::endl;
	       throw(std::invalid_argument(errMsg));
	    }
	    for (size_t i = 0; i < fsize; ++i){
	       *(fbegin+i) /= valarray[i];
	    }
	    return (*this);
	 }

   private:
      //Basic constuctor is private
      ArraySlice () {}
      const size_t fsize; //The size should not change
      T * fbegin;
      const T * fbeginConst;
};

template<typename T>
std::ostream & operator << (std::ostream &os, const ArraySlice<T> & arr){
   os<<"[";
   for (size_t i = 0; i < arr.size()-1; ++i){
      os<<arr[i]<<",";
   }
   os<<arr[arr.size()-1]<<"]";
   return os;
}
#endif
