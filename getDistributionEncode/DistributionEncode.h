#include <armadillo>
#include <map>
#include <vector>
#include <list>
#include <set>

//#include "mlpack/methods/sparse_coding/random_initializer.hpp"
//#include "mlpack/methods/sparse_coding/sparse_coding.hpp"
#define  shiftN 12
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706
#define Bin  0
#define Cosine 1
//using namespace mlpack::sparse_coding;
using namespace arma;
using namespace std;

typedef struct clusterEntry
{
   int Rank;
   std::vector <int> IDArray;
   std::vector <int> PatternIDs;
}clusterEntry;



class DistributionEncode{

     public:
     std::vector<clusterEntry> clusterResult;

     double tolerance =0.01;
     int MaxDepth = 1024;
     bool rankSkip = true;
     bool IsRead = false;
     int minsupport=1;
     double RBFThreshold = 0.1;
     double Diffusion=32;
     unsigned int HistType = Bin;
     void Vectorization(sp_mat & ConnectionMatrix, int N, std::vector< std::vector< std::complex<double> > > &DistributionEncodeList);
     void Vectorization_bins(sp_mat & ConnectionMatrix, int N, std::vector< std::vector< std::complex<double> > > &DistributionEncodeList);
     void TableReadWithName(FILE *ip, std::vector< std::vector < unsigned int > >&Sample,std::vector < std::string >&SampleName);
     void GenoType(FILE * Input,FILE * Output);
     void PairwiseGenotype(FILE * Input,FILE * Output);
     void PairwiseEncode(mat&InputM,FILE * Output,int K);
     void RBFEncode( mat&InputM,FILE * Output,int K,int CompsN );
     void Encode(mat&InputM,FILE * Output,int K );
     void MeanEncode(mat &InputM,FILE * Output);

};





