# DistributionEncode

# Dependencies
  Armadillo package and mlpack package are necessary. 
# Compile
     make
  The execute file will be generated.

# Input file format 
  The input file format is csv. A delimited text file that uses a comma to separate values without any head name and row name.
 
  Each group contain many entries(such as cells), each entry assign an entry ID. Each entry has many attribute values(such as gene expression).
 
 The first column is the group ID, the second column is entry ID, remaining columns are the atrribute values(such as genes expressions in each cell).
 
 group ID  entry ID  atrribute value1  atrribute value2   ...  atrribute valueN ...
 
    1,        1,               20,            16,        ...   6,   ...
    1,        2,               19,            21,        ...   101, ...
    1,        3,               21,            50,        ...   12,  ...
    1,        4,               20,            20,        ...   10, ...
    1,        5,               18,            7,         ...   48, ...
    1,        6,               46,            12,        ...   29, ...
    1,        7,               13,             7,        ...   35, ...
    1,        8,                4,            26,        ...   0, ...
    1,        9,                7,            19,        ...   7, ...
    1,        10,              35,            24,        ...   72, ...
    1,        11,              26,            10,        ...   1, ...
    1,        12,              58,             1,        ...   29, ...
    1,        13,              15,            20,        ...   31, ...
    1,        14,               22,            4,        ...   0, ...
    1,        15,               52,           22,        ...   5, ...
    1,        16,               37,           18,        ...   0, ...
    2,        17,               27,           66,        ...   35, ... 
    2,        18,               41,           65,        ...  14, ... 
    2,        19,               12,           37,        ... 33, ...
    2,        20,               32,           20,        ... 42, ...
    2,        21,               34,           23,        ... 24, ...
    2,        22,               25,           32,        ... 23, ...
    2,        23,                8,           36,        ... 69, ...
    .........
   # Comands
   
      ./DistributionEncode Encode InputFile Outputfile 
      
  InputFile is a csv file described in the "Input the file format" section. The output is a table with the number of columns N and the number of rows M. N is the number of orthognal basis functions, M = n\*m n is the number of samples and m is the number of attributes (such as genes).  Each line of this table is a vector represent the entries distribution for an attribute of a sample. n lines are as a group with have the same attribute.
  
  
   
    ./DistributionEncode RBFEncode InputFile Outputfile [ComponentsNumber]
    
  InputFile is a csv file described in the "Input the file format" section. ComponentsNumber is an integer indicating the number of consider components. If this parameter is specified, the output is a table with the number of columns N and the number of rows M.
N is the number of consider components. M = n\*m n is the number of samples and m is the number of attributes (such as genes) in each sample. Each line of this table is the vector of components represent the distribution of entries for an attribute of a sample. n lines are as a group with have the same attribute.     
 If ComponentsNumber parameter is not specified, the output is not a regular table, each line has different columns(the length of vector is variable). But lines belong to the same attribute have the same column.
 
    ./DistributionEncode RBFEncode -s sampleN InputFile Outputfile [ComponentsNumber]
 
 InputFile is a csv file as the output file of  "DistributionEncode Encode InputFile Outputfile" command，sampleN is the number of samples. ComponentsNumber is an integer indicating the number of consider components. If this parameter is specified, the output is a table with the number of columns N and the number of rows M. N is the number of consider components. M = n\*m n is the number of samples and m is the number of attributes (such as genes) in each sample. Each line of this table is the vector of components represent the distribution of entries for an attribute of a sample. n lines are as a group with have the same attribute.     
 If ComponentsNumber parameter is not specified, the output is not a regular table, each line has different columns(the length of vector is variable). But lines belong to the same attribute have the same column.
 
 
    ./DistributionEncode MeanEncode InputFile Outputfile  
   
  InputFile is a csv file described in the "Input the file format" section. The output is a table with the number of columns N and the number of rows M. N is the number of samples and the M is the number of attributes. In this table each element is the average value of entries belong to same attribute and the same sample. 
   
   
   
    ./DistributionEncode PairwiseEncode InputFile Outputfile   
     
   InputFile is a csv file described in the "Input the file format" section. The output is a table with the number of columns N and the number of rows M. N is the number of pairs of samples (n\*(n-1)/2 n the number of samples). M is the number of attributes. In this table ,each element is the difference of entries distribution between two samples.    
   
    ./DistributionEncode PairwiseEncode -s sampleN InputFile Outputfile
   
   InputFile is a csv file as the output file of  "DistributionEncode Encode InputFile Outputfile" command，sampleN is the number of samples. The output is a table with the number of columns N and the number of rows M. N is the number of pairs of samples (n\*(n-1)/2 n the number of samples). M is the number of attributes. In this table ,each element is the difference of entries distribution between two samples.  
    
   
    ./DistributionEncode  PairwiseGenotype InputFile Outputfile
   
  inputfile 
      
      s1  s2  s3 
      AA  TT  AT
      .. 
  Outputfile 
  
      s1_s2  s1_s3   s2_s3  
       2       1  　    1
      .. 
   The Output File is correspond to Output File of Pairwise Encode command. 
      
      
# config file 
      
      Tolerance = 0.01
      VectorSize =   512
      RBFThreshold = 0.1
      
  ”VectorSize“  is the size of the encoding vector, ”tolerance“ is the resolution of orthogonal basis function. RBF Threshold is used to select the components, if the "ComponentsNumber" parameter is not specified in "RBFEncode" command,  the rato of eigenvalue for selected component and the sum of all eigenvalues should be above this value. 

    -config   configPATH/config 
    
 If the config file and execute file are not locate in the same file directory, the path should be specified.  
     
