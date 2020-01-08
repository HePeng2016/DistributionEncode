#include "DistributionEncode.h"




void DistributionEncode::Vectorization(sp_mat & ConnectionMatrix, int N, std::vector< std::vector< std::complex<double> > > &DistributionEncodeList)
{



    std::vector <unsigned int> GroupCount;
    DistributionEncodeList.resize(ConnectionMatrix.n_rows);
    GroupCount.resize(ConnectionMatrix.n_rows);
    int PN;
    double GaussianSum =0;
    double Energy =0;
    std::complex<double> Image(0,1);
    double t_0 = 0.01;
    double t_step=3.1415926*2/(N);// sqrt(N*tolerance);
    std::vector <double> Sum;





    for(int i=0;i<(ConnectionMatrix.n_rows);i++)
    {
        DistributionEncodeList[i].resize(N);
        GroupCount[i]=0;

        for(int j =0;j<N;j++)
        {
          DistributionEncodeList[i][j]=0.0;
        }

    }


    for(PN=0;PN<N;PN++)
    {
      if( exp(-pow((PN*Diffusion/N),2)/2)<0.00000001)
      break;
    }

    for(int i=0;i<N;i++)
    {

    }
    for(int i=0;i<clusterResult.size();i++)
    {




        if(clusterResult[i].IDArray.size()<2)
          continue;


        Energy = ConnectionMatrix(clusterResult[i].IDArray[0],clusterResult[i].IDArray[1]);
        for( int n=0;n<PN;n++)
        {
           DistributionEncodeList[clusterResult[i].IDArray[0]][n] = DistributionEncodeList[clusterResult[i].IDArray[0]][n]+ exp((-pow((n*Diffusion/N),2)/2))*exp((round(Energy/tolerance)+Diffusion)*Image*(t_step*n));
        }
        GroupCount[clusterResult[i].IDArray[0]]= GroupCount[clusterResult[i].IDArray[0]]+1;
     }

      std::complex<double> GaussianBias;



     for(int i=0;i<(ConnectionMatrix.n_rows);i++)
     {
         GaussianBias = 0.0;
         for( int n=0;n<N;n++)
         {
            GaussianBias = GaussianBias+DistributionEncodeList[i][n];
            DistributionEncodeList[i][n] = DistributionEncodeList[i][n]/(((double)GroupCount[i])+0.0000000000000001);
         }

            GaussianBias = GaussianBias/(((double)GroupCount[i])+0.0000000000000001);
            DistributionEncodeList[i][0]= DistributionEncodeList[i][0]-GaussianBias;

     }


}


void DistributionEncode::TableReadWithName(FILE *ip, std::vector< std::vector < unsigned int > >&Sample,std::vector < std::string >&SampleName)
{
  std::vector<double>Temp;
  char buffer[1024];
  unsigned int entry = 1;
  int N =0;
  int M =0;

   std::string tempword;

    do{


	    fscanf(ip,"%[\b|\t]*",buffer);
	    fscanf(ip,"%s",&buffer);
            N++;
            fscanf(ip,"%[\b|\t]*",buffer);
	} while(!fscanf(ip,"%[\n]",buffer));

     rewind(ip);
     while(!feof(ip))
	{
		  entry = fgetc(ip);
		  if(entry=='\n')
		  {
		      M++;
		  }
	 }
   rewind(ip);
   Sample.resize(M-1);
   for(int i=0;i<(M-1);i++)
   {
     Sample[i].resize(N);
   }
   SampleName.resize(N);


    for(int j=0;j<N;j++)
    {
              	fscanf(ip,"{\b|\t|\n}*",buffer);
			    fscanf(ip,"%s",buffer);
			    SampleName[j].assign(buffer);
    }
    for(int i=0;i<M-1;i++)
    {	 for(int j=0;j<N;j++)
		 {
			    double value;
              	fscanf(ip,"{\b|\t|\n}*",buffer);
			    fscanf(ip,"%s",buffer);
				Sample[i][j]=(unsigned int)buffer[0];
				Sample[i][j]=(Sample[i][j]<<shiftN)|(unsigned int)buffer[1];
		 }
	 }
}


void DistributionEncode::GenoType(FILE * Input,FILE * Output)
  {


        std::vector< std::vector < unsigned int > >Sample;
        std::vector < std::string >SampleName;
        int  Vec[512];

         TableReadWithName(Input,Sample,SampleName);

         int I = 0;
         int SampleID_Back=0;

        for(int i=0;i<SampleName.size();i++)
        {
          fprintf(Output,"%s",SampleName[i].c_str());
          fprintf(Output,"\t");
        }
       fprintf(Output,"\n");

         for(int n=0;n<Sample.size();n++)
         {

            Vec['A']=0;
            Vec['T']=0;
            Vec['C']=0;
            Vec['G']=0;

            for(int i=0;i<Sample[n].size();i++)
            {
                int A_i = Sample[n][i]&0xFFF;
                int B_i = (Sample[n][i]>>shiftN)&0xFFF;
                Vec[A_i]=Vec[A_i]+1;
                Vec[B_i]=Vec[B_i]+1;
            }

            int Max=0;
            short int MaxI;
            if(Vec['A']>Max)
            {
                MaxI = 'A';
                Max  = Vec['A'];
            }

            if(Vec['T']>Max)
            {
                MaxI = 'T';
                Max  = Vec['T'];
            }

            if(Vec['C']>Max)
            {
                MaxI = 'C';
                Max  = Vec['C'];
            }
            if(Vec['G']>Max)
            {
                MaxI = 'G';
                Max  = Vec['G'];
            }

            for(int i=0;i<Sample[n].size();i++)
            {
                double score;

                int A_i = Sample[n][i]&0xFFF;
                int B_i = (Sample[n][i]>>shiftN)&0xFFF;
                if( (A_i==MaxI)&&(B_i==MaxI))
                {
                    score=0.0;

                }else if( (A_i==MaxI)||(B_i==MaxI) )
                {
                    score=1.0;
                }else
                 {
                    score=2.0;
                 }
                fprintf(Output,"%lf",score);
                fprintf(Output,"\t");

            }
                fprintf(Output,"\n");
           }

    }

void DistributionEncode::PairwiseGenotype(FILE * Input,FILE * Output)
{


        std::vector< std::vector < unsigned int > >Sample;
        std::vector < std::string >SampleName;


        TableReadWithName(Input,Sample,SampleName);


        for(int i=0;i<SampleName.size();i++)
        {
            for(int j=0;j<i;j++)
            {
                fprintf(Output,"%s_",SampleName[i].c_str());
                fprintf(Output,"%s\t",SampleName[j].c_str());
            }

        }
        fprintf(Output,"\n");

         for(int n=0;n<Sample.size();n++)
          {
            for(int i=0;i<Sample[n].size();i++)
              for(int j=0;j<i;j++)
               {
                   double score;
                   int A_i = Sample[n][i]&0xFFF;
                   int B_i = (Sample[n][i]>>shiftN)&0xFFF;
                   int A_j = Sample[n][j]&0xFFF;
                   int B_j = (Sample[n][j]>>shiftN)&0xFFF;


                     if(Sample[n][i]==Sample[n][j])
                     {
                        score=0.0;

                     }else if ((A_i==A_j)||( B_i==B_j ))
                     {
                        score=1.0;

                     }else
                     {
                        score=2.0;
                     }
                      fprintf(Output,"%lf\t",score);
               }
                      fprintf(Output,"\n");
           }
    }

void DistributionEncode::PairwiseEncode(mat&InputM,FILE * Output,int K)
{


  if(IsRead ==false)
  {
        double SampleID = InputM(InputM.n_rows-1,0);
        double CellID = InputM(InputM.n_rows-1,1);
        sp_mat EnergyMatrix;
        mat OutputM;
        cx_vec TempY;
        TempY.set_size(K);

        EnergyMatrix.resize((SampleID),(CellID));
        clusterResult.resize(InputM.n_rows);

        for(int i=0;i<clusterResult.size();i++)
        {
            clusterResult[i].IDArray.resize(2);
            clusterResult[i].IDArray[1]=(InputM(i,1)-1);
            clusterResult[i].IDArray[0]=(InputM(i,0)-1);
        }

        OutputM.resize(EnergyMatrix.n_rows,K);


        for(int n=2;n<InputM.n_cols;n++)
        {
            for(int i=0;i<InputM.n_rows;i++)
           {
              EnergyMatrix(( InputM(i,0)-1),(InputM(i,1)-1))=InputM(i,n);
           }
           std::vector< std::vector< std::complex<double> > > DistributionEncodeList;
           Vectorization(EnergyMatrix,K,DistributionEncodeList);


            for(int  i=0;i<DistributionEncodeList.size();i++)
           {
             for(int j =0;j<DistributionEncodeList[i].size();j++)
             {
                 TempY[j]=DistributionEncodeList[i][j];
             }
                 TempY=fft(TempY);

                for(int j =0;j<DistributionEncodeList[i].size()-1;j++)
                {
                   double value = abs(real(TempY[j])/(K/2.0));
                   if(value>0.000000001)
                   OutputM(i,j)= sqrt(value);
                   else
                   OutputM(i,j)= 0;
                }
           }





            for(int  i=0;i<DistributionEncodeList.size();i++)
            {
                for(int j =0;j<i;j++)
                {
                  double x=norm((OutputM.row(j) - OutputM.row(i)),2);
                  fprintf(Output,"%lf\t",x);
                }
            }

                 fprintf(Output,"\n");

          }
  }else
  {
         for(int n=0;n<(InputM.n_rows-K+1);)
         {
            for(int  i=n;i<n+K;i++)
            {
                for(int j =n;j<i;j++)
                {
                  double x=norm((InputM.row(j) - InputM.row(i)),2);
                  fprintf(Output,"%lf\t",x);
                }
            }
            fprintf(Output,"\n");
            n=n+K;
         }

  }
}

void DistributionEncode::RBFEncode( mat&InputM,FILE * Output,int K,int CompsN )
{

   if(IsRead ==false)
   {
        mat OutputM;
        double SampleID = InputM(InputM.n_rows-1,0);
        double CellID = InputM(InputM.n_rows-1,1);
        cx_vec TempY;
        TempY.set_size(K);

        sp_mat EnergyMatrix;
        EnergyMatrix.resize((SampleID),(CellID));
        clusterResult.resize(InputM.n_rows);

        for(int i=0;i<clusterResult.size();i++)
        {
            clusterResult[i].IDArray.resize(2);
            clusterResult[i].IDArray[1]=(InputM(i,1)-1);
            clusterResult[i].IDArray[0]=(InputM(i,0)-1);
        }

        OutputM.resize(EnergyMatrix.n_rows,K);


         vec eigval;
         mat eigvec;
         mat KernalM;
         double Comps =3;
         eigvec.resize(EnergyMatrix.n_rows,Comps);
         KernalM.resize(EnergyMatrix.n_rows,EnergyMatrix.n_rows);

         for( int i=0;i<EnergyMatrix.n_rows;i++)
         {
            KernalM(i,i)=1.0;
         }


        for(int n=2;n<InputM.n_cols;n++)
        {
            for(int i=0;i<InputM.n_rows;i++)
           {
              EnergyMatrix(( InputM(i,0)-1),(InputM(i,1)-1))=InputM(i,n);
           }
           std::vector< std::vector< std::complex<double> > > DistributionEncodeList;

           Vectorization(EnergyMatrix,K,DistributionEncodeList);



           for(int  i=0;i<DistributionEncodeList.size();i++)
           {
                for(int j =0;j<DistributionEncodeList[i].size();j++)
                {
                   TempY[j]=DistributionEncodeList[i][j];
                }
                   TempY=fft(TempY);

                for(int j =0;j<DistributionEncodeList[i].size()-1;j++)
                {
                   double value = abs(real(TempY[j])/(K/2.0));
                   if(value>0.000000001)
                   OutputM(i,j)= sqrt(value);
                   else
                   OutputM(i,j)= 0;
                }
           }





            for(int  i=0;i<DistributionEncodeList.size();i++)
            {
                for(int j =0;j<i;j++)
                {
                  double x=sum(square(OutputM.row(j) - OutputM.row(i)));
                  x=exp(x*(-0.5));
                  KernalM(i,j)=x;
                  KernalM(j,i)=x;
                }
            }

               eig_sym(eigval,eigvec,KernalM);

               if(CompsN>0)
               {
                 Comps=CompsN;
               }else
               {

                   for(int i =(KernalM.n_rows-1);i>=0;i--)
                   {

                     if(eigval(i)/KernalM.n_rows < RBFThreshold)
                     {
                        Comps = KernalM.n_rows -i;
                        break;
                     }
                   }

               }

               for(int i=0;i<KernalM.n_rows;i++)
               {

                 for(int j=0;j<Comps;j++)
                 {
                   fprintf(Output,"%lf	",eigvec(i,(KernalM.n_rows-(j+1))));
                 }
                 fprintf(Output,"\n");
               }
          }
    }else
    {

        vec eigval;
        mat eigvec;
        mat KernalM;
        double Comps =3;
        eigvec.resize(K,Comps);
        KernalM.resize(K,K);

         for( int i=0;i<K;i++)
         {
            KernalM(i,i)=1.0;
         }

         for(int i=0;i<InputM.n_rows;i++)
            for(int j=0;j<InputM.n_cols;j++)
         {
             double value = abs(InputM(i,j));

              if(value>0.000000001)
              InputM(i,j)= sqrt(value);
              else
              InputM(i,j)= 0;
         }



        for(int n=0;n<(InputM.n_rows-K+1);)
        {

            for(int  i=n;i<n+K;i++)
            {
                for(int j =n;j<i;j++)
                {
                  double x=sum(square(InputM.row(j) - InputM.row(i)));
                  x=exp(x*(-0.5));
                  KernalM(i-n,j-n)=x;
                  KernalM(j-n,i-n)=x;
                }
            }

               eig_sym(eigval,eigvec,KernalM);

               if(CompsN>0)
               {
                 Comps=CompsN;
               }else
               {

                   for(int i =(KernalM.n_rows-1);i>=0;i--)
                   {

                     if(eigval(i)/KernalM.n_rows < RBFThreshold)
                     {
                        Comps = KernalM.n_rows -i;
                        break;
                     }
                   }

               }

               for(int i=0;i<KernalM.n_rows;i++)
               {

                 for(int j=0;j<Comps;j++)
                 {
                   fprintf(Output,"%lf	",eigvec(i,(KernalM.n_rows-(j+1))));
                 }
                 fprintf(Output,"\n");
               }
            n=n+K;
        }

    }

}

void DistributionEncode::Encode(mat&InputM,FILE * Output,int K )
{


        double SampleID = InputM(InputM.n_rows-1,0);
        double CellID = InputM(InputM.n_rows-1,1);


        sp_mat AttributesMatrix;
        AttributesMatrix.resize((SampleID),(CellID));
        cx_vec TempY;
        TempY.set_size(K);


        for(int n=2;n<InputM.n_cols;n++)
       {

           for(int i=0;i<InputM.n_rows;i++)
           {
             AttributesMatrix(( InputM(i,0)-1),(InputM(i,1)-1))=InputM(i,n);
           }

             clusterResult.resize(InputM.n_rows);

           for(int i=0;i<clusterResult.size();i++)
           {
             clusterResult[i].IDArray.resize(2);
             clusterResult[i].IDArray[1]=(InputM(i,1)-1);
             clusterResult[i].IDArray[0]=(InputM(i,0)-1);
           }


           std::vector< std::vector< std::complex<double> > > DistributionEncodeList;

           Vectorization(AttributesMatrix,K,DistributionEncodeList);



           for(int  i=0;i<DistributionEncodeList.size();i++)
           {


             for(int j =0;j<DistributionEncodeList[i].size();j++)
             {
                 TempY[j]=DistributionEncodeList[i][j];
             }
                TempY=fft(TempY);

           for(int j =0;j<DistributionEncodeList[i].size()-1;j++)
           {
               fprintf(Output,"%lf,",real(TempY[j])/(K/2.0));
           }
               fprintf(Output,"%lf  ",real(TempY[DistributionEncodeList[i].size()-1])/(K/2.0));
               fprintf(Output,"\n");

           }

      }

    }


void DistributionEncode::MeanEncode(mat &InputM,FILE * Output)
{

        double SimpleID = InputM(InputM.n_rows-1,0);


        mat SummaryMatrix;
        SummaryMatrix.zeros((SimpleID),InputM.n_cols-2);
        mat Count;
        Count.zeros((SimpleID),1);





        for(int i=0;i<InputM.n_rows;i++)
        {
           for(int n=2;n<InputM.n_cols;n++)
          {
            SummaryMatrix((InputM(i,0)-1),(n-2))=InputM(i,n)+SummaryMatrix((InputM(i,0)-1),(n-2));
          }

           Count((InputM(i,0)-1),0)=Count((InputM(i,0)-1),0)+1;
        }


        for(int i=0;i<SummaryMatrix.n_cols;i++)
        {
           for(int j=0;j<SummaryMatrix.n_rows;j++)
           {
              SummaryMatrix(j,i)=SummaryMatrix(j,i)/Count(j,0);

           }
        }

        for(int i=0;i<SummaryMatrix.n_cols;i++)
        {
            for(int j=0;j<SummaryMatrix.n_rows;j++)
            {
                fprintf(Output,"%lf",SummaryMatrix(j,i));
                fprintf(Output,"\t");
            }
                fprintf(Output,"\n");
         }

  }
