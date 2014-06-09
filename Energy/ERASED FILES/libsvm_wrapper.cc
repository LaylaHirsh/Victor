#include "libsvm_wrapper.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <Debug.h>
using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

char *line;
int max_line_len = 1024;
//double lower=-1.0,upper=1.0,y_lower,y_upper;
double y_lower,y_upper;
int y_scaling = 0;
//double *feature_max;
//double *feature_min;
double y_max = -1e10;
double y_min = 1e10;
int max_index;
int ifc; 
/************************************************************/

// Constructor. Set parameters to default values
LibSVMWrapper::LibSVMWrapper() {
  prob.x = NULL;
  prob.y = NULL; 
  x_space = NULL;
  model = NULL;

  param.svm_type = EPSILON_SVR;
  param.kernel_type = RBF; 
  param.degree = 3;
  param.gamma = 0;	// 1/k
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = 1;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;
  //cross_validation = 0;
}

// Constructor. Read model from file
LibSVMWrapper::LibSVMWrapper(const string& modelname, 
			     const string& restore_filename):
  lower(-1.0), upper(1.0) {

  prob.x = NULL;
  prob.y = NULL;
  x_space = NULL;
  if((model = svm_load_model(modelname.c_str())) == 0) {
    fprintf(stderr, "can't open model file %s\n", modelname.c_str());
    exit(1);
  }

  // read file for scaling features
  FILE *fp = fopen(restore_filename.c_str(),"r");
  int idx, c;
  double fmin, fmax;
  
  if(fp==NULL)
    ERROR("Can't open file.", exception);

  if((c = fgetc(fp)) == 'y') {
	  ifc=fscanf(fp, "%lf %lf\n", &y_lower, &y_upper);
	  ifc=fscanf(fp, "%lf %lf\n", &y_min, &y_max);
    y_scaling = 1;
  } else
    ungetc(c, fp);
  
  if (fgetc(fp) == 'x') {
    ifc=fscanf(fp, "%lf %lf\n", &lower, &upper);
    while(fscanf(fp,"%d %lf %lf\n",&idx,&fmin,&fmax)==3) {
      feature_min.push_back(fmin);
      feature_max.push_back(fmax);
      //if(idx<=max_index) {
      //feature_min[idx] = fmin;
      //feature_max[idx] = fmax;
      //}
    }
  }
  fclose(fp);
}

LibSVMWrapper::~LibSVMWrapper() {
  svm_destroy_model(model);
  svm_destroy_param(&param);

  free(prob.y);
  free(prob.x);
  free(x_space);
}

struct svm_parameter LibSVMWrapper::getParameter() { return param; }

void LibSVMWrapper::setParameter(struct svm_parameter parameter) { 
  param = parameter; 
}

void LibSVMWrapper::load(const string& modelname) {
  if((model = svm_load_model(modelname.c_str())) == 0) {
    fprintf(stderr, "can't open model file %s %d\n", modelname.c_str(), ifc);
    exit(1);
  }
}

void LibSVMWrapper::save(const string& modelname) {
  svm_save_model(modelname.c_str(), model);
}

// TODO
// Get and set external training parameters (e.g. C, ...)
// For the moment those parameters are set externally
// and then train in launched.
void LibSVMWrapper::train(const string& datasetfile) {
  if(prob.y != NULL) free(prob.y);
  if(prob.x != NULL) free(prob.x);
  if(x_space != NULL) free(x_space);

  read_problem(datasetfile.c_str());
  const char* error_msg = svm_check_parameter(&prob, &param);  
  if(error_msg) {
    fprintf(stderr,"Error: %s\n",error_msg);
    exit(1);
  }

  model = svm_train(&prob, &param);
}

// void LibSVMWrapper::test() {}

// TODO
double LibSVMWrapper::compute_margin(struct svm_node*) { return .0; }

double LibSVMWrapper::operator()(const vector<double>& features) {
  vector<double>::size_type inputdim = features.size();
  assert(inputdim == feature_min.size());

  struct svm_node *x = (struct svm_node *) malloc((inputdim+1)*sizeof(struct svm_node));
  // build (sparse?) representation of the input
  for(vector<double>::size_type j=0; j<inputdim; ++j) {
    x[j].index = j+1;
    x[j].value = output(j, features[j]); 
    //cout << '(' << x[j].value << ',' << features[j] << ") " << flush;
  }
  //  cout << endl;
  x[inputdim].index = -1;
  
  double tmp = svm_predict(model, x);
  if (tmp > 1.0)
    tmp = 1.0;
  else if (tmp < 0.0)
    tmp = 0.0;
  return tmp;
}

/***************** private methods *******************/

// read in a problem (in svmlight format)
void LibSVMWrapper::read_problem(const char *filename) {
  int elements, max_index, i, j;
  FILE *fp = fopen(filename,"r");
	
  if(fp == NULL) {
    fprintf(stderr,"can't open input file %s\n",filename);
    exit(1);
  }

  prob.l = 0;
  elements = 0;
  while(1) {
    int c = fgetc(fp);
    switch(c) {
    case '\n':
      ++prob.l;
      // fall through,
      // count the '-1' element
    case ':':
      ++elements;
      break;
    case EOF:
      goto out;
    default:
      ;
    }
  }
 out:
  rewind(fp);

  prob.y = Malloc(double,prob.l);
  prob.x = Malloc(struct svm_node *,prob.l);
  x_space = Malloc(struct svm_node,elements);

  max_index = 0;
  j=0;
  for(i=0;i<prob.l;i++) {
    double label;
    prob.x[i] = &x_space[j];
    ifc=fscanf(fp,"%lf",&label);
    prob.y[i] = label;
    while(1) {
      int c;
      do {
	c = getc(fp);
	if(c=='\n') goto out2;
      } while(isspace(c));
      ungetc(c,fp);
      ifc=fscanf(fp,"%d:%lf",&(x_space[j].index),&(x_space[j].value));
      ++j;
    }	
  out2:
    if(j>=1 && x_space[j-1].index > max_index)
      max_index = x_space[j-1].index;
    x_space[j++].index = -1;
  }

  if(param.gamma == 0)
    param.gamma = 1.0/max_index;
  
  fclose(fp);
}

// output scaled feature
double LibSVMWrapper::output(int index, double value) {
  /* skip single-valued attribute */
  if(feature_max[index] == feature_min[index])
    return value;

  if(value == feature_min[index])
    value = lower;
  else if(value == feature_max[index])
    value = upper;
  else
    value = lower + (upper-lower) * 
      (value-feature_min[index])/
      (feature_max[index]-feature_min[index]);
  
  return value;

  //if(value != 0)
  //printf("%d:%g ",index, value);
}
