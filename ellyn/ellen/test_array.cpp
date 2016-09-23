#include <iostream>
#include <vector>
using namespace std;

void test_array(float **features, int N, int D){

  std::vector< std::vector< float > > X(N,vector<float>());

  for (unsigned int i = 0; i< X.size(); ++i){
    X[i].assign(features[i],features[i]+D);
    std::cout << "features[i]:" << *(features[i]) << "\n";
    std::cout << "features[i]+D:" << *(features[i]+D) << "\n";
    std::cout << "x[" << i << "]";
    for (auto j : X[i])
      std::cout << j << ",";
    std::cout << "\n";
  }
  std::cout << "\nsize of X:" << X.size() << "x" << X[0].size() << "\n";

}

int main(){
  float arr[2][3] = {{0.4, 5.2, 5.19}, {6.54, 3.21, 0.10}};
  vector<float *> b;

  for (size_t i = 0; i < 2; ++i)
    b.push_back(arr[i]);

  test_array(b.data(),2,3);

}
