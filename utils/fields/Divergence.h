#ifndef __DIVERGENCE_H__
#define __DIVERGENCE_H__

#include "ScalarFieldExpressions.h"
#include "VectorField.h"
#include <cstddef>
#include <functional>

namespace fdaPDE{
namespace core{

  // a functor representing the divergence of a VectorField<N>
  template <int N>
  class Divergence : public FieldExpr<Divergence<N>>{
  private:
    // the function encoding the divergence operator
    std::function<double(SVector<N>)> div_;

  public:
    Divergence(const VectorField<N>& field) {
      // set up the functor encoding the divergence operator
      std::function<double(SVector<N>)> div = [field](SVector<N> x) -> double {
	double result = 0;
	for(std::size_t i = 0; i < N; ++i){
	  result += field[i].derive()(x)[i]; // sum of derivatives evaluated at point x
	}
	return result;
      };
      div_ = div;
    };

    double operator()(const SVector<N>& x) const { return div_(x); };
  };

  template <int N>
  Divergence<N> div(const VectorField<N>& field){
    return Divergence<N>(field);
  }
  
}};
#endif // __DIVERGENCE_H__
