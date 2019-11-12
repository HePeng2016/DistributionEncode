# DistributionEncode

# Dependencies
  Armadillo package and mlpack package are necessary. 

# Input file format 
  The input file format is csv. A delimited text file that uses a comma to separate values without any head name and row name.
 
  Each group contain many entries(such as cells), each entry assign an entry ID. Each entry has many attribute values(such as gene expression).
 
 The first column is the group ID, the second column is entry ID, remaining columns are the atrribute values(genes expressions in each cell).
 
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
    
