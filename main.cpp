#include  "getDistributionEncode/DistributionEncode.h"
using namespace std;


int main( int argc , char *argv[] )
{

     DistributionEncode  Test;
     Test.tolerance=0.01;
     double K =1024;
     FILE * configip=NULL;
     char buffer[1024];


    for(int i =0;i<argc-1;i++)
    {
       if( strcmp(argv[i],"-config")==0)
       {
           configip = fopen(argv[i+1],"r");
       }
    }

   if( configip == NULL)
    {
        configip = fopen("config","r");
    }

   if(configip != NULL )
   {
     while (!feof(configip))
    {
        fscanf(configip,"%[\b|\t]*",buffer);
        fscanf(configip,"%s",buffer);

        if(strcmp(buffer,"Tolerance")==0)
        {
            fscanf(configip,"%[\b|\t]*",buffer);
            fscanf(configip,"%s",buffer);

            if(strcmp(buffer,"=")==0)
            {
               fscanf(configip,"%lf",&Test.tolerance);
            }
               fscanf(configip,"%[\b|\t|\n]",buffer);
         continue;
       }

       if(strcmp(buffer,"VectorSize")==0)
       {
            fscanf(configip,"%[\b|\t]*",buffer);
            fscanf(configip,"%s",buffer);

            if(strcmp(buffer,"=")==0)
            {
               fscanf(configip,"%lf",&K);
            }
               fscanf(configip,"%[\b|\t|\n]",buffer);
         continue;
       }

       if(strcmp(buffer,"RBFThreshold")==0)
       {
            fscanf(configip,"%[\b|\t]*",buffer);
            fscanf(configip,"%s",buffer);

            if(strcmp(buffer,"=")==0)
            {
               fscanf(configip,"%lf",&Test.RBFThreshold);
            }
               fscanf(configip,"%[\b|\t|\n]",buffer);
         continue;
       }


    }
   }


if(argc<2)
{
  printf("GenoType\n");
  printf("MeanEncode\n");
  printf("Encode\n");
  printf("PairwiseGenotype\n");
  printf("PairwiseEncode\n");
  printf("RBFEncode\n");

}

for(int i =0;i<argc;i++)
{


  if ( strcmp(argv[i],"GenoType") ==0 )
  {

        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("GenoType InputFile Outputfile \n");
           return 0;
        }


        FILE * ip1 = fopen(argv[i+1],"r");
        FILE * ip2 = fopen(argv[i+2],"w+");
        if( ip1==NULL )
        {
           printf("The input file does not exist.\n");
           return 0;
        }
        if( ip2==NULL )
        {
           printf("The output File cannot be created.\n");
           return 0;
        }

        Test.GenoType(ip1,ip2);
        fclose(ip2);
        fclose(ip1);
        return 0;
    }


  if ( strcmp(argv[i],"MeanEncode") ==0 )
  {

        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("MeanEncode InputFile Outputfile \n");
           return 0;
        }

        mat InputM;
        if( InputM.load(argv[i+1],csv_ascii)==false)
        {
           printf("The Input File cannot be loaded as a matrix.\n");
           return 0;
        }
        FILE * ip2 = fopen(argv[i+2],"w+");
        if( ip2==NULL )
        {
           printf("The output File cannot be created.\n");
           return 0;
        }

        Test.MeanEncode(InputM,ip2);
        fclose(ip2);
        return 0;

  }


  if ( strcmp(argv[i],"PairwiseGenotype") ==0 )
  {

        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("PairwiseGenotype InputFile Outputfile \n");
           return 0;
        }

        FILE * ip1 = fopen(argv[i+1],"r");
        FILE * ip2 = fopen(argv[i+2],"w+");


        if( ip1==NULL )
        {
           printf("The input file does not exist.\n");
           return 0;
        }
        if( ip2==NULL )
        {
           printf("The output File cannot be created.\n");
           return 0;
        }


        Test.PairwiseGenotype(ip1,ip2);
        fclose(ip2);
        fclose(ip1);
        return 0;
    }


  if ( strcmp(argv[i],"PairwiseEncode") ==0 )
  {

        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("PairwiseGenotype InputFile Outputfile \n");
           printf("PairwiseGenotype -m InputFile Outputfile \n");
           return 0;
        }

        mat InputM;

        if(InputM.load(argv[i+1],csv_ascii)==false)
        {
           printf("The Input File cannot be loaded as a matrix.\n");
           return 0;
        }

       FILE * ip2 = fopen(argv[i+2],"w+");

       if( ip2==NULL )
       {
           printf("The output File cannot be created.\n");
           return 0;
       }
       Test.PairwiseEncode(InputM,ip2,K);
       fclose(ip2);
       return 0;
  }

  if ( strcmp(argv[i],"RBFEncode") ==0 )
  {
        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("RBFEncode InputFile Outputfile [ComponentsNumber]\n");
           printf("RBFEncode -m InputFile Outputfile [ComponentsNumber]\n");
           return 0;
        }

        mat InputM;
        if( InputM.load(argv[i+1],csv_ascii) ==false)
        {
           printf("The Input File cannot be loaded as a matrix.\n");
           return 0;
        }


        FILE * ip2 = fopen(argv[i+2],"w+");

        if( ip2==NULL )
        {
           printf("The output File cannot be created.\n");
           return 0;
        }

        int compsN;
        if( argv[i+3]!= NULL)
         compsN = atoi(argv[i+3]);
         else
         compsN = -1;
        Test.RBFEncode(InputM,ip2,K,compsN);
        fclose(ip2);
        return 0;
    }



  if ( strcmp(argv[i],"Encode") ==0 )
  {

        if(argv[i+1]==NULL||argv[i+2]==NULL)
        {
           printf("Encode InputFile Outputfile\n");
           return 0;
        }

        mat InputM;

        if( InputM.load(argv[i+1],csv_ascii) ==false)
        {
           printf("The Input File cannot be loaded as a matrix.\n");
           return 0;
        }
        InputM.load(argv[i+1],csv_ascii);
        FILE * ip3 = fopen(argv[i+2],"w+");
        if( ip3==NULL )
        {
           printf("The output File cannot be created.\n");
           return 0;
        }

        Test.Encode(InputM,ip3,K);
        fclose(ip3);
        return 0;
    }


}



}



