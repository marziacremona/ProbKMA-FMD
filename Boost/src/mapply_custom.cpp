#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
// [[Rcpp::depends(RcppParallel)]]
struct MapplyCustom : public Worker {
  const Function fun;
  const List args;
  const bool simplify;
  const bool use_names;
  RObject result;
  List results;
  
  MapplyCustom(Function fun, List args, bool simplify, bool use_names)
    : fun(fun), args(args), simplify(simplify), use_names(use_names) 
  {
    results = List(args.size());
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      RObject res = fun(args[i]);
      results[i] = res;
    }
  }
  
  void join(const MapplyCustom& rhs) 
  {
    for (std::size_t i = 0; i < results.size(); ++i) 
    {
      RObject lhs_result = results[i];
      RObject rhs_result = rhs.results[i];
      
      if (simplify) 
      {
        if (use_names) 
        {
          CharacterVector lhs_names = lhs_result.attr("names");
          CharacterVector rhs_names = rhs_result.attr("names");
          int total_length = lhs_names.size() + rhs_names.size();
          CharacterVector combined_names(total_length);
          for (int j = 0; j < lhs_names.size(); ++j) {
            combined_names[j] = lhs_names[j];
          }
          for (int j = 0; j < rhs_names.size(); ++j) {
            combined_names[lhs_names.size() + j] = rhs_names[j];
          }
          List list_result = List::create(lhs_result, rhs_result);
          list_result.attr("names") = combined_names;
          results[i] = list_result;
        } else 
        {
          results[i] = List::create(lhs_result, rhs_result);
        }
      }else 
      {
        results[i] = List::create(lhs_result, rhs_result);
      }
    }
  }
};

// [[Rcpp::export]]
List mapply_custom(Function fun, List args, bool simplify = true, bool use_names = true) 
{
  MapplyCustom mapplyFun(fun, args, simplify, use_names);
  parallelFor(0, args.size(), mapplyFun); //parallelize a loop across multiple threads or processors
  
  return mapplyFun.results;
}
