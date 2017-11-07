#ifndef _parameteric_curves_utils_file_io_hpp
#define _parameteric_curves_utils_file_io_hpp

#include <fstream>  /// to read text file
#include <Eigen/Dense>

#define MAXBUFSIZE  ((int) 1000000)

namespace parametriccurves {
namespace utils {

const Eigen::MatrixXd readMatrixFromFile(const std::string& filename)
{
  int cols = 0, rows = 0;
  double buff[MAXBUFSIZE];
  
  // Read numbers from file into buffer.
  std::ifstream infile;
  infile.open(filename.c_str());
  while (! infile.eof())
  {
    std::string line;
    getline(infile, line);
    //          std::cout<<"Read line "<<rows<<"\n";
    
    int temp_cols = 0;
    std::stringstream stream(line);
    while(! stream.eof())
      stream >> buff[cols*rows+temp_cols++];
    
    if (temp_cols == 0)
      continue;
    
    if (cols == 0)
      cols = temp_cols;
    else if(temp_cols!=cols && !infile.eof())
    {
      std::cout<<"Error while reading matrix from file, line "<<rows<<" has "<<temp_cols<<" columnds, while preceding lines had "<<cols<<" columnds\n";
      std::cout<<line<<"\n";
      break;
    }
    
    rows++;
    if((rows+1)*cols>=MAXBUFSIZE)
    {
      std::cout<<"Max buffer size exceeded ("<<rows<<" rows, "<<cols<<" cols)\n";
      break;
    }
  }
  infile.close();
  rows--;
  
  // Populate matrix with numbers.
  Eigen::MatrixXd result(rows,cols);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      result(i,j) = buff[ cols*i+j ];
  return result;
}

} //namespace utils
} //namespace parametriccurves

#endif //_CLASS_EXACTCUBIC
