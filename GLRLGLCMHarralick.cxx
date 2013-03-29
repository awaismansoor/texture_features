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

#include <itkScalarImageTextureCalculator.h>
#include "itkMinimumMaximumImageCalculator.h"


#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
/*
The code calculates:
11 GLRL FEATURES
   Short Run Emphasis (SRE)
   Long Run Emphasis (LRE)
   Gray-Level Nonuniformity (GLN)
   Run Length Nonuniformity (RLN)
   Run Percentage (RP)
   Low Gray-Level Run Emphasis (LGRE)
   High Gray-Level Run Emphasis (HGRE)
   Short Run Low Gray-Level Emphasis (SRLGE)
   Short Run High Gray-Level Emphasis (SRHGE)
   Long Run Low Gray-Level Emphasis (LRLGE)
   Long Run High Gray-Level Emphasis (LRHGE) 
   
6 GLCM HARALICK FEATURES
	Energy
	Entropy
	InverseDifferenceMoment
	Inertia
	ClusterShade
	ClusterProminence
   
8 MOMEMTS
	  Mean
      Mean Dev
      Std Dev
      Variance
      Skewness
      Kurtosis
      Min
      Max

	  
 ---------------------------------------------
 Author:
 ---------------------------------------------
    (C)Awais Mansoor <awais.mansoor@gmail.com>
    

*/

using namespace cv;
typedef itk::Image< unsigned char, 3 > ImageType;



void HaralickFeatures( const ImageType* image, double stats[6])
{
		typedef itk::Statistics::ScalarImageTextureCalculator<ImageType> GLCMGeneratorType;
		GLCMGeneratorType::Pointer texture_measure = GLCMGeneratorType::New();
		texture_measure->FastCalculationsOn();
		texture_measure->SetNumberOfBinsPerAxis( 16 );
		texture_measure->SetInsidePixelValue( 255 );
					
		texture_measure->SetInput( image );
		texture_measure->Compute();
		const GLCMGeneratorType::FeatureValueVector* output = texture_measure->GetFeatureMeans();
		for(unsigned int ii = 0; ii < output->size(); ++ii)		
			stats[ii]=(*output)[ii];

						
					
		std::cout << "Energy = "<<(*output)[0]<< std::endl;
		std::cout << "Entropy = "<<(*output)[1]<< std::endl;
		std::cout << "InverseDifferenceMoment = "<<(*output)[2]<< std::endl;
		std::cout << "Inertia = "<<(*output)[3]<< std::endl;
		std::cout << "ClusterShade = "<<(*output)[4]<< std::endl;
		std::cout << "ClusterProminence = "<<(*output)[5]<< std::endl;

}

void momentosEstatOpenCV(IplImage *src, double stats[8]) {

	/*
      stats[0] = mean
      stats[1] = mean dev
      stats[2] = std dev
      stats[3] = variance
      stats[4] = skewness
      stats[5] = kurtosis
      stats[6] = min
      stats[7] = max
  */
  printf("momentosEstatOpenCV\n");

  double 	epR = 0.0, epG = 0.0, epB = 0.0, sumR, sumG, sumB, pR, pB, pG;
  int i, j, w, h, N;
	CvScalar pixel;

  
  	for (i = 0; i < 8; i++) 
			stats[i] = 0.0;
	stats[6] = 255.0;
  

  // Step 1: Compute mean, min and max
  sumR = 0.0;
  sumG = 0.0;
  sumB = 0.0;
  w = src->width;
  h = src->height;
  N = w*h;
  for (i = 0; i < h; i++) 
    for (j = 0; j < w; j++) {
			pixel = cvGet2D(src, i, j);
			sumB += pixel.val[0];

      stats[6] = (stats[6] < pixel.val[0]) ? stats[6] : pixel.val[0];
      
      stats[7] = (stats[7] > pixel.val[0]) ? stats[7] : pixel.val[0];
      
    }
  stats[0] = (double)sumB/N;


  // Step 2: compute second, third and fourth moments
  //         from mean deviation
  for (int i = 0; i < h; i++) 
    for (int j = 0; j < w; j++) {
			pixel = cvGet2D(src, i, j);
			sumB = pixel.val[0] - stats[0];

			stats[1] += fabs(sumB);

      epB += sumB;
      epR += sumR;
      epG += sumG;
			stats[3] += (pB = sumB*sumB);
			stats[4] += (pB *= sumB);

			stats[5] += (pB *= sumB);

    }
  stats[1] /= N;
  
  stats[3] = (stats[3] - epB*epB/N)/(N-1);
  
  stats[2] = sqrt(stats[3]);
  
  if (stats[3] > 0) {
     stats[4] /= (N*(stats[3])*(stats[2]));
     stats[5] = (stats[5])/(N*stats[3]*stats[3]) - 3.0;
  } else
     printf("No skew/kurtosis. Variance = 0.\n");

    std::cout<<"Mean = "<<stats[0]<<std::endl;
	std::cout<<"Mean Dev. = "<<stats[1]<<std::endl;
	std::cout<<"Std. Dev. = "<<stats[2]<<std::endl;
	std::cout<<"Variance = "<<stats[3]<<std::endl;
	std::cout<<"Skewness = "<<stats[4]<<std::endl;
	std::cout<<"Kurtosis = "<<stats[5]<<std::endl;
	std::cout<<"Min. = "<<stats[6]<<std::endl;
	std::cout<<"Max. ="<<stats[7]<<std::endl;

} 


void momentosEstatITK( const ImageType* image, double stats[8]) {
  /*
      stats[0] = mean
      stats[1] = mean dev
      stats[2] = std dev
      stats[3] = variance
      stats[4] = skewness
      stats[5] = kurtosis
      stats[6] = min
      stats[7] = max
  */

printf("momentosEstatITK\n");


  double 	ep = 0.0, sum, p;
  int i, j, w, h, N;


  
  	for (i = 0; i < 8; i++) 
			stats[i] = 0.0;
	stats[6] = 255.0;
  

  // Step 1: Compute mean, min and max
  sum = 0.0;
  w=image->GetRequestedRegion().GetSize(0);
  h=image->GetRequestedRegion().GetSize(1);

  ImageType::IndexType pixelIndex;
  
  N = w*h;
  double val=0.0;
  for (i = 0; i < h; i++) 
    for (j = 0; j < w; j++) {			
			
			pixelIndex[0]=j;
			pixelIndex[1]=i;
			pixelIndex[2]=0;

			sum+=image->GetPixel(pixelIndex);

      stats[6] = (stats[6] < image->GetPixel(pixelIndex)) ? stats[6] : image->GetPixel(pixelIndex);
      
      stats[7] = (stats[7] > image->GetPixel(pixelIndex)) ? stats[7] : image->GetPixel(pixelIndex);
      
    }

  stats[0] = (double)sum/N;

  


  // Step 2: compute second, third and fourth moments
  //         from mean deviation
  for (int i = 0; i < h; i++) 
    for (int j = 0; j < w; j++) {

			pixelIndex[0]=j;
			pixelIndex[1]=i;
			pixelIndex[2]=0;

			sum = image->GetPixel(pixelIndex) - stats[0];

			stats[1] += fabs(sum);

			ep += sum;

			stats[3] += (p = sum*sum);
			stats[4] += (p *= sum);

			stats[5] += (p *= sum);

    }


  stats[1] /= N;
  
  stats[3] = (stats[3] - ep*ep/N)/(N-1);
  
  stats[2] = sqrt(stats[3]);
  
  if (stats[3] > 0) {
     stats[4] /= (N*(stats[3])*(stats[2]));
     stats[5] = (stats[5])/(N*stats[3]*stats[3]) - 3.0;
  } else
     printf("no skew/kurtosis. Variance = 0.");

    std::cout<<"Mean = "<<stats[0]<<std::endl;
	std::cout<<"Mean Dev. = "<<stats[1]<<std::endl;
	std::cout<<"Std. Dev. = "<<stats[2]<<std::endl;
	std::cout<<"Variance = "<<stats[3]<<std::endl;
	std::cout<<"Skewness = "<<stats[4]<<std::endl;
	std::cout<<"Kurtosis = "<<stats[5]<<std::endl;
	std::cout<<"Min. = "<<stats[6]<<std::endl;
	std::cout<<"Max. ="<<stats[7]<<std::endl;


}
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

  
  double Stats[25];
  double stats2[8];

  Stats[0]=SRE;
  Stats[1]=LRE;
  Stats[2]=GLN;
  Stats[3]=RLN;
  Stats[4]=RP;
  Stats[5]=LGRE;
  Stats[6]=HGRE;
  Stats[7]=SGLGE;
  Stats[8]=SRHGE;
  Stats[9]=LRLGE;
  Stats[10]=LRHGE;

  
  IplImage *imagem = 0;
  imagem = cvLoadImage(argv[1], 1);
  momentosEstatITK( reader->GetOutput(), stats2);

  for(i=0;i<8; i++)
	  Stats[11+i]=stats2[i];



						
		

  //momentosEstatOpenCV(imagem, stats2);
  double stats3[6];
  HaralickFeatures( reader->GetOutput(), stats3);

  for(i=0;i<6; i++)
	  Stats[19+i]=stats3[i];

  for(i=0;i<25; i++)
	  std::cout<<Stats[i]<<std::endl;
 


  return EXIT_SUCCESS;
}


