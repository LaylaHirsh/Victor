#ifndef _LIBSVM_WRAPPER_H_
#define _LIBSVM_WRAPPER_H_

#include "svm.h"
#include <string>
#include <vector>
using namespace std;

class LibSVMWrapper {
  struct svm_parameter param;
  struct svm_problem prob; // set by read_problem
  struct svm_model *model;
  struct svm_node *x_space;

  double lower, upper;
  vector<double> feature_min;
  vector<double> feature_max;
  
  void read_problem(const char *filename);
  double output(int, double);

 public:
  LibSVMWrapper();
  LibSVMWrapper(const std::string&, const std::string&);
  ~LibSVMWrapper();

  struct svm_parameter getParameter();
  void setParameter(struct svm_parameter);

  void train(const std::string&);
  //void test();
  double compute_margin(struct svm_node*); // TODO
  double operator()(const std::vector<double>&);
  void load(const std::string&);
  void save(const std::string&);

};

#endif // _LIBSVM_WRAPPER_H_
