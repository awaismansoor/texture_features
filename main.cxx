#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkAddImageFilter.h"
#include <iostream>
#include "math.h"
#include "string.h"
#include "malloc.h" 

#include <itkScalarImageToGreyLevelCooccurrenceMatrixGenerator.h> 
#include "itkMinimumMaximumImageCalculator.h"


#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//The program 

using namespace cv;


int RoundtoNearestInt(double input)
{
return (floor(input+0.5));
}

Mat mesh_grid1(Mat c_vector, Mat r_vector)
{
	transpose(r_vector, r_vector);
	
	Mat MeshGrid=Mat::zeros(r_vector.size[0], c_vector.size[1], CV_32F);

	for(int i=0; i<r_vector.size[0]; i++)
		{
			
			c_vector.copyTo(MeshGrid.row(i));
			//std::cout<<MeshGrid.row(i)<<std::endl;
		}

	return MeshGrid;

}

Mat mesh_grid2(Mat c_vector, Mat r_vector)
{
	transpose(r_vector, r_vector);
	
	Mat MeshGrid=Mat::zeros(r_vector.size[0], c_vector.size[1], CV_32F);

	for(int j=0; j<c_vector.size[1]; j++)
		{
			
			r_vector.copyTo(MeshGrid.col(j));
			//std::cout<<MeshGrid.row(i)<<std::endl;
		}

	return MeshGrid;

}

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile thresholdLow thresholdHigh" << std::endl;
      return EXIT_FAILURE;
    }
  
  //ITK settings
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef unsigned char OutPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutPixelType, Dimension > OutImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
    
  //Filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
   
  //Parameters
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  float thresholdLow = atof(argv[3]);
  float thresholdHigh = atof(argv[4]);

  //Pipeline
  try
    {
      reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  size = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = size[0];
  pCol = size[1];
  pSli = size[2]; 
  ImageType::RegionType region;
  region.SetSize( size );
  //Allocate new image
  ImageType::Pointer image = ImageType::New();
  image->SetRegions( region );
  image->SetSpacing( spacing );
  image->SetOrigin( origin );
  image->SetDirection( direction );
  image->Allocate();  

  typedef itk::MinimumMaximumImageCalculator <ImageType>
          ImageCalculatorFilterType;

  ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
  imageCalculatorFilter->SetImage(reader->GetOutput());
  imageCalculatorFilter->Compute();

  

  int NL=8;
  double slope = (NL-1)/(imageCalculatorFilter->GetMaximum()-imageCalculatorFilter->GetMinimum());
  double intercept =  1 - (slope*(imageCalculatorFilter->GetMinimum()));


  typedef itk::MultiplyByConstantImageFilter<ImageType, double, ImageType> MultiplyImageFilterType;
  MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(reader->GetOutput());
  multiplyImageFilter->SetConstant(slope);

  multiplyImageFilter->Update();

  typedef itk::AddConstantToImageFilter <ImageType, double, ImageType> AddImageFilterType;
  AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput(multiplyImageFilter->GetOutput());  
  addImageFilter->SetConstant(intercept);
  addImageFilter->Update();

  
  //Read image
  ImageType::IndexType pixelIndex;
  int i, j, k, value;
  for(i=0;i<pRow;i++)
    for(j=0;j<pCol;j++)      
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=0;
	  value = RoundtoNearestInt(addImageFilter->GetOutput()->GetPixel(pixelIndex));
	  image->SetPixel(pixelIndex, value);
	  

	}

ImageType::IndexType pixelIndex2;

int* index=(int*) malloc (512*sizeof(int));
int* len=(int*) malloc (512*sizeof(int));
int* val=(int*) malloc (512*sizeof(int));



Mat oneglrlm = Mat::zeros(NL, pCol, CV_32F);



int ii=0;

 for(i=0; i<pRow; i++)
 {
	ii=0;
	Mat temp = Mat::zeros(NL, pCol, CV_32F);

	for(j=0;j<pCol;j++)      
	{
	  pixelIndex[0]=i;
	  pixelIndex[1]=j;
	  pixelIndex[2]=0;

	  pixelIndex2[0]=i;
	  pixelIndex2[1]=j+1;
	  pixelIndex2[2]=0;
	  
	  if (image->GetPixel(pixelIndex)!=image->GetPixel(pixelIndex2))
	  {
		  index[ii]=j;
		  val[ii]=image->GetPixel(pixelIndex);		  
		  ii++;
	  }

	  


	}
	 
	 len[0]=index[0]+1;


	  for(k=1; k<ii; k++)
	  {
		  len[k]=index[k]-index[k-1];		  
		  len[ii]=pCol-1-index[k];		  
	  
	  }

	  for(k=0; k<ii; k++)
	  {		  		  
		  temp.row(val[k]-1).col(len[k]-1)=temp.row(val[k]-1).col(len[k]-1)+1;
	  }

	 oneglrlm+= temp;

	  //std::cout<<oneglrlm.size()<<std::endl;
	
 }

 
 int numStats = 11;
 Mat stats = Mat::zeros(1, numStats, CV_32F);

 Mat c_vector=Mat::zeros(1, oneglrlm.size[0], CV_32F);
 Mat r_vector=Mat::zeros(1, oneglrlm.size[1], CV_32F);
 
 Mat p_g=Mat::ones(1, oneglrlm.size[0], CV_32F);
 Mat p_r=Mat::ones(oneglrlm.size[1], 1, CV_32F);
 
 for(i=0; i<oneglrlm.size[0]; i++)
	 c_vector.col(i)=i+1;
 
 for(j=0; j<oneglrlm.size[1]; j++)
	 r_vector.col(j)=j+1;


 p_g=p_g*oneglrlm;
 p_r=oneglrlm*p_r;




 // Total number of runs
 int N_runs = sum(p_g)[0];

 //total number of elements
 int N_tGLRLM = oneglrlm.size[0]*oneglrlm.size[1];

 transpose(p_r,p_r);
 Mat SRE_vector=Mat::zeros(1, oneglrlm.size[0], CV_32F);
 multiply(c_vector, c_vector, SRE_vector);
 divide(p_r, SRE_vector, SRE_vector);
 double SRE=sum(SRE_vector)[0]/N_runs;

 Mat LRE_vector=Mat::zeros(1, oneglrlm.size[0], CV_32F);
 multiply(c_vector, c_vector, LRE_vector);
 multiply(p_r, LRE_vector, LRE_vector);
 double LRE=sum(LRE_vector)[0]/N_runs;

 Mat GLN_vector=Mat::zeros(1, oneglrlm.size[1], CV_32F);
 multiply(p_g, p_g, GLN_vector); 
 double GLN=sum(GLN_vector)[0]/N_runs;

 Mat RLN_vector=Mat::zeros(1, oneglrlm.size[0], CV_32F);
 multiply(p_r, p_r, RLN_vector); 
 double RLN=sum(RLN_vector)[0]/N_runs;

 double RP=(double)N_runs/N_tGLRLM;

 Mat LGRE_vector=Mat::zeros(1, oneglrlm.size[1], CV_32F);
 multiply(r_vector, r_vector, LGRE_vector);
 divide(p_g, LGRE_vector, LGRE_vector);
 double LGRE=sum(LGRE_vector)[0]/N_runs;

 Mat HGRE_vector=Mat::zeros(1, oneglrlm.size[1], CV_32F);
 multiply(r_vector, r_vector, HGRE_vector);
 multiply(p_g, HGRE_vector, HGRE_vector);
 double HGRE=sum(HGRE_vector)[0]/N_runs;

 Mat c_matrix=mesh_grid1(c_vector, r_vector);
 Mat r_matrix=mesh_grid2(c_vector, r_vector);

 Mat SGLGE_matrix=Mat::zeros(oneglrlm.size[1], oneglrlm.size[1], CV_32F);
 multiply(r_matrix, c_matrix, SGLGE_matrix);
 multiply(SGLGE_matrix, SGLGE_matrix, SGLGE_matrix);
 transpose(SGLGE_matrix, SGLGE_matrix);
 divide(oneglrlm, SGLGE_matrix, SGLGE_matrix);
 Mat temp=Mat::ones(SGLGE_matrix.size[1], 1, CV_32F);
 double SGLGE=sum(SGLGE_matrix*temp)[0]/N_runs;

 Mat SRHGE_matrix=Mat::zeros(oneglrlm.size[1], oneglrlm.size[1], CV_32F);
 multiply(r_matrix, r_matrix, SRHGE_matrix);
 transpose(oneglrlm, oneglrlm);
 multiply(oneglrlm, SRHGE_matrix, SRHGE_matrix); 
 divide(SRHGE_matrix, c_matrix, SRHGE_matrix);
 divide(SRHGE_matrix, c_matrix, SRHGE_matrix);
 temp=Mat::ones(SRHGE_matrix.size[1], 1, CV_32F);
 double SRHGE=sum(SRHGE_matrix*temp)[0]/N_runs;

 Mat LRLGE_matrix=Mat::zeros(oneglrlm.size[1], oneglrlm.size[1], CV_32F);
 multiply(c_matrix, c_matrix, LRLGE_matrix);
 multiply(oneglrlm, LRLGE_matrix, LRLGE_matrix); 
 divide(LRLGE_matrix, r_matrix, LRLGE_matrix);
 divide(LRLGE_matrix, r_matrix, LRLGE_matrix);
 temp=Mat::ones(LRLGE_matrix.size[1], 1, CV_32F);
 double LRLGE=sum(LRLGE_matrix*temp)[0]/N_runs;


 Mat LRHGE_matrix=Mat::zeros(oneglrlm.size[1], oneglrlm.size[1], CV_32F);
 multiply(c_matrix, c_matrix, LRHGE_matrix);
 multiply(oneglrlm, LRHGE_matrix, LRHGE_matrix);
 multiply(LRHGE_matrix, r_matrix, LRHGE_matrix);
 multiply(LRHGE_matrix, r_matrix, LRHGE_matrix);
 temp=Mat::ones(LRHGE_matrix.size[1], 1, CV_32F);
 double LRHGE=sum(LRHGE_matrix*temp)[0]/N_runs; 

  	
  std::cout << "SRE = "<< SRE<<std::endl;
  std::cout << "LRE = "<< LRE<<std::endl;
  std::cout << "GLN = "<< GLN<<std::endl;
  std::cout << "RLN = "<< RLN<<std::endl;
  std::cout << "RP = "<< RP<<std::endl;
  std::cout << "LGRE = "<< LGRE<<std::endl;
  std::cout << "HGRE = "<< HGRE<<std::endl;
  std::cout << "SGLGE = "<< SGLGE<<std::endl;
  std::cout << "SRHGE = "<< SRHGE<<std::endl;
  std::cout << "LRLGE = "<< LRLGE<<std::endl;
  std::cout << "LRHGE = "<< LRHGE<<std::endl;

  

  return EXIT_SUCCESS;
}


