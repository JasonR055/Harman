/**
 * \file Harman/src/Harman.cpp 
 * \date 15/01/2015
 * \author   Josh Bowden, CSIRO
 * \details Harman eResearch Collaboration project ERRFP-263 (https://jira.csiro.au/browse/ERRFP-263)
 * Based on Matlab code developed by Yalchin Oytam
 * Stash site: https://<userid><at>stash.csiro.au/scm/~bow355/harman_r.git
 * 
 * Install from source:
 * Windows:
 *  Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source Harman
 * Linux:
 *   R CMD INSTALL --preclean --build Harman -l ~/R/library
 * 
 * DOxygen : BUILTIN_STL_SUPPORT    = YES
 * Use devtools::use_rcpp() to create some C++ friendly package details
 * 
 * The C++ functionality has doxygen based markup.
 * The doxygen configuration file (docs/Harman.config) uses the 'dot' tool for
 * interactive UML type graphs of code structure, which may only be availble through Linux
 * versions and can be disabled through setting HAVE_DOT = NO in the config file.
 *  To make the C++ documentation:
 *    <from bash prompt>
 *      cd <git repo>/harman_r/docs
 *      doxygen Harman.config
 * 
 * To create the R documentation via roxygen:
 *  <from R prompt>
 *    library(roxygen2)
 *    roxygen2::roxygenise()

 Windows:
 library(Harman)
 batches <- read.table("F:\\projects\\harman_git\\data\\harman\\Rob_Dunne_info.csv", header=TRUE, sep=",", row.names="Sample")
 expt_<- batches$Treatment
 batch_ <- as.factor(batches$Batch)
 datamatrix <- read.table("F:\\projects\\harman_git\\data\\harman\\Rob_Dunne_data_original.csv", header=TRUE, sep=",", row.names="probeID")
 system.time(harmanresults  <- harman(as.matrix(datamatrix), expt=expt_, batch=batch_, strict=FALSE, printInfo=T) )
 
 library(Harman)
 data(NPM)
 expt <- (npm.info$Treatment)
 batch <- (npm.info$Batch)
 system.time(harmanresults <- harman(as.matrix(npm.data), expt=expt, batch=batch, limit=0.95, printInfo = TRUE))
 system.time(harmanresults <- harman(as.matrix(npm.data), expt=expt, batch=batch, limit=0.95, printInfo = TRUE, forceRand=T))

 
 Linux:
 batches <- read.table("/home/bow355/ACP_HARMAN2015/harman_git/data/Rob_Dunne_info.csv", header=TRUE, sep=",", row.names="Sample")
 datamatrix <- read.table("/home/bow355/ACP_HARMAN2015/harman_git/data/Rob_Dunne_data_original.csv", header=TRUE, sep=",", row.names="probeID")
 expt_<- batches$Treatment
 batch_ <- batches$Batch
 system.time(harmanresults  <- harman(as.matrix(datamatrix), expt=expt_, batch=batch_) )

 
 
*/
 
#define _USE_MATH_DEFINES
// Standard template library headers
#include <set>
#include <map>
#include <list>
#include <vector>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <iostream>
#include <string>
#include <iomanip>

// only need for debugging?
#include <fstream>
#include <sstream>

#ifdef _OPENMP
//#ifdef SUPPORT_OPENMP
    #include <omp.h>
#endif

//// [[Rcpp::depends()]]
//#include <progress.hpp>

//static Progress *  progmon = NULL ;

// for linux timing function clock_gettime()  Requires -lrt for the "real time" library
#include <time.h>  
#if _LINUX
    #include <unistd.h> // check that _POSIX_TIMERS is > 1
#endif 

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif


// #define _ARRAYALIGNEMT (4)
#ifdef _WIN32
#define _DATAPATH "F:\\projects\\harman_git\\data\\harman\\"
#else
#define _DATAPATH "/home/bow355/ACP_HARMAN2015/harman_git/data/"
#endif


#ifdef _WIN32
#define strtoll     _strtoi64
#define strtoull    _strtoui64
#if defined (_MSC_VER) 
#if defined (WIN32)
// http://stackoverflow.com/questions/2170385/c-math-functions round() not defined in VS2010 math library
double round(double val)
{    
    return floor(val + 0.5);
}
#endif
#endif
#endif



#if (_USE_RCPP==1)
//#define _RANDFUNCTION_ (unif_rand() * ((unsigned int)-1))
#define _PRINTERROR Rcpp::Rcerr
#define _PRINTSTD Rcpp::Rcout
#else
//#define _RANDFUNCTION_ rand()
#define _PRINTERROR std::cerr
#define _PRINTSTD std::cout
#endif


#if !defined (_MSC_VER)
// include this if we want to use __float128 and its special sprintf version
// also need to add -lquadmath to PKG_LIBS// 34 significant digits
// __float128 gives 34 significant digits
// N.B. nchoosek(300,130) ~ 10^84 combinations (order does not matter, however no repetition of a choice is allowed)
// may need to look into using GNU Multiple Prec. Library: https://gmplib.org/
/*#if !defined (__INTEL_COMPILER)
extern "C" {
#include <quadmath.h>
}
#endif
*/


// Rcpp R interface library
#include <Rcpp.h>



// This is the main function for the Harman package functionality
RcppExport SEXP HarmanMain(SEXP pc_scores_rIN, SEXP group_rIN, SEXP limit_rIN, SEXP numrepeats_rIN, SEXP randseed_rIN, SEXP forceRandom_rIN, SEXP printInfoIN) ;
  
// These are all testing functions
//#define _PRINTSTD Rcpp::Rcout 
//RcppExport SEXP HarmanTestDataStructures(SEXP scores_r, SEXP batch_r, SEXP limit_r) ;
//RcppExport SEXP HarmanTestBatchTreatmentExtraction(SEXP scores_r, SEXP batch_r, SEXP limit_, SEXP treatmentnum_r, SEXP batchnum_r ) ;
//RcppExport SEXP test_nchoosek_rcpp_factorials_class(SEXP N_IN, SEXP K_IN) ;
//RcppExport SEXP TestHarmanCSimulateBatchDistribution(SEXP scores_r, SEXP batch_r, SEXP limit_r) ;
//RcppExport SEXP test_randomSelctFromIntegerVect(SEXP inputDataIN, SEXP numselectIN, SEXP numrepsIN, SEXP seedIN) ;
//#if !defined (__INTEL_COMPILER)
//RcppExport SEXP test_randomSelctBiasVsUnbias(SEXP numValsIN, SEXP startRangeIN, SEXP endRangeIN, SEXP seedIN) ;
#endif 



#if defined (_MSC_VER)
// This library is available for windows but not through Microsoft
// Library build is available here:		F:\aa_ASC_projects\ACP_epiGPU_onf\getopt_mb_uni_vc10\getopt_mb_uni_vc10\getopt_mb_uni_vc10.sln
// getops.h header is available here:					F:\aa_ASC_projects\ACP_epiGPU_onf\getopt_mb_uni_vc10\getopt_mb_uni_vc10\getopt_mb_uni_vc10_dll
// getops.dll is available here: 	F:\aa_ASC_projects\ACP_epiGPU_onf\getopt_mb_uni_vc10\getopt_mb_uni_vc10\bin\release\dll
// with License: LGPL
#undef _UNICODE
#include <getopt.h>
#endif

// #define _TESTRANDOMVERSION 1
#include "CMapSelectKFromN.h"
#include "CSelectRandom.h"
int TOTAL_ITERATIONS ;



// using namespace Rcpp; // i'll be explicit for every namespace being used
 
/*!  \brief Used to specify the layout of an experiments data matrix - [samples x variables] of [variables x samples] 
   * \details SamplesAsRows is a recognised default for a data matrix - each row is a sample with variables in columns.
   * SamplesAsRows also means that on a C/C++ row major system that all variables for a particular sample are adjacent (linear) in memory
   *  
   * Example data layout:
   * S = sample, V = variable
   * 
   * SamplesAsRows:
   *        V1  V2  V3  V4
   *    S1  -   -   -   -
   *    S2  -   -   -   - 
   *    S3  -   -   -   -
   * 
   * SamplesAsCols:
   *        S1  S2  S3 
   *    V1  -   -   - 
   *    V2  -   -   - 
   *    V3  -   -   - 
   *    V4  -   -   -
   */
enum DataOrientation {SamplesAsRows=0, SamplesAsCols=1} ;


/*! \brief Inputs and takes a copy of the complete scores data and alows access to scores of a defined PC 
   */
class  CPCAScoresArray
{
private:
  
public:
  std::vector<double> _scores ;   /*! this will act as a 2D array storage  */
  size_t _numPCS ;                /*! the number of columns in _scores */
  size_t _numSamples ;            /*! the number of rows in _scores */
  DataOrientation _dataOriented ; /*! if (_dataOriented == ) rows = _numSamples, cols = _numPCS */

  /*! **Creator method**
   * \details Copies input data over so that scores for a particular PC are available in a row 
   */
  CPCAScoresArray(std::vector<double> scoresIN, size_t samplesIN, size_t PCsIN,  DataOrientation dataOriented)
  {

      _scores.reserve(scoresIN.size()) ;
      _scores.assign(&scoresIN[0], &scoresIN[0]+scoresIN.size());  //TODO: Test this assignement method
      this->_numPCS = PCsIN ;
      this->_numSamples = samplesIN ;
      this->_dataOriented = dataOriented ;
  }
  
  /*! 
   * \details Selects the scores data of a particular PC and returns it in a std::vector<double>
	 * \param PCWantedIN is a 0 based index into the stored _scores data
	 * \returns a std::vector<double> of PC scores
   */
  std::vector<double> * GetPCData(size_t PCWantedIN) // PCWantedIN is 0 based. returns a std::vector<double> of PC scores
  {
    std::vector<double> * retVect = (std::vector<double> *) new std::vector<double>(this->_numSamples) ; // constructor - create the return vector the size of the number of samples present
    double * ptr_retVect =  retVect->data() ;
    
    if (this->_dataOriented == SamplesAsRows)
    {
      size_t start = PCWantedIN ;
      for (size_t i = 0 ; i < this->_numSamples; i++)
      {
        ptr_retVect[i] = this->_scores[start+(i*this->_numPCS)] ; // move through data 'column', in steps of numPCs
      }
    }
    else
    {
			size_t start = PCWantedIN * this->_numSamples ;
      retVect->assign(&_scores[start], &_scores[start]+this->_numSamples);      
    }
    
    return (retVect) ;
  } 
  
  /*! 
   * \details Sets the scores data of a particular PC 
   * \param   PCToSetIN is a 0 based index into the stored _scores data
   * \param   scoreDataIN is the vector of scores data used to replace the PCToSetIN scores
   * \returns void
   */
  void SetPCData(size_t PCToSetIN, std::vector<double> * scoreDataIN ) // PCToSetIN is 0 based. returns a std::vector<double> of PC scores
  {
    double * ptr_scoreDataIN =  scoreDataIN->data() ;
    
    if (this->_dataOriented == SamplesAsRows)
    {
      size_t start = PCToSetIN ;
      for (size_t i = 0 ; i < this->_numSamples; i++)
      {
         this->_scores[start+(i*this->_numPCS)] = ptr_scoreDataIN[i] ; // move through data 'column', in steps of numPCs
      }
    }
    else
    {
      size_t start = PCToSetIN * this->_numSamples ;
      for (size_t i = 0 ; i < this->_numSamples; i++)
      {
        this->_scores[start+i] = ptr_scoreDataIN[i] ; // move through data 'row', in steps of sizeof(double)
      }
    }
    
  } 
  
  /*! 
   * \details Sets the scores data of a particular PC 
   * \param   PCToSetIN is a 0 based index into the stored _scores data
   * \param   scoreDataIN is the vector of scores data used to replace the PCToSetIN scores
   * \returns void
   */
  void SetPCData_TEST(size_t PCToSetIN, std::vector<double> * scoreDataIN ) // PCToSetIN is 0 based. returns a std::vector<double> of PC scores
  {
    //double * ptr_scoreDataIN =  scoreDataIN->data() ;
    
    if (this->_dataOriented == SamplesAsRows)
    {
     // _PRINTSTD << "Testing : SamplesAsRows " << std::endl ;
      size_t start = PCToSetIN ;
      for (size_t i = 0 ; i < this->_numSamples; i++)
      {
        this->_scores[start+(i*this->_numPCS)] = (double)PCToSetIN ; // move through data 'column', in steps of numPCs
      }
    }
    else
    {
     // _PRINTSTD << "Testing : SamplesAsRows " << std::endl ;
      size_t start = PCToSetIN * this->_numSamples ;
      for (size_t i = 0 ; i < this->_numSamples; i++)
      {
        this->_scores[start+i] = (double) PCToSetIN  ; // move through data 'row', in steps of sizeof(double)
      }
    }
    
  } 
  
  /*! \brief Returns a pointer to the first element of double vector
   */
  double * ReturnPointerToVectorData( void )
  {
    return &(this->_scores[0]) ;
  }
  
};



/*! \brief Holds the scores (T) data for rapid access to specific batch and treatment samples
 Class to hold the scores (T) data of a single PC as a list of list of pointers to vectors of doubles used for rapid access to specific batch and treatment samples.
* A CExperimentStructure class is needed to create the _llpvd_T structure from an input vector of scores data
*/
class  CExperimentData
{
  std::list< std::list<std::vector<double>* > > _llpvd_T ; /*! _llpvd_T is a list of lists of pointers to vectors of doubles (llpvd)
                                                          Inner list contains vectors of values (for example PCA scores) for a particular treatment type,
                                                          with a seperate vector containing scores for that treatment present in each batch
                                                          and the outer list is a list of each treatment type.
																													\verbatim                                                      
	Example data layout:
	t = treatment b = batch
	list[ list( t1_b1<2.1,3.2>,						t1_b2<3.2,0.8>,			t1_b3<0.1,3.1> ),      
				list( t2_b1<1.8,1.8>,						t2_b2<2.1,-3.5>,		t2_b3<-2.2,1.2> ),   
				list( t3_b1<2.5,3.4,2.1,-3.5>,	t3_b2<5.1,-2.2>,		t3_b3<3.2,-0.8> ),
				list( t4_b1<>,									t4_b2<2.1,3.2>,			t4_b3<3.2,-2.8> ) ]
																													 \endverbatim                             
		\see Return_T_WithVariableData() 
																													*/  

public:

	CExperimentData(size_t num_trIN, size_t num_baIN) 
	{
		this->Initialise_T_VectorsOfSamples( num_trIN,  num_baIN) ;
	}

	~CExperimentData( void ) 
	{
		this->Free_T_VectorsOfVariableData() ;
	}


	std::list< std::list<std::vector<double>* > > * ReturnData( void ) 
	{
		return &_llpvd_T ;
	}

  /*!
   * \details  Allocates the vectors<doubles> of T_llvd structure that can be filled with samples of a particular treatment and batch combination
   */
   void Initialise_T_VectorsOfSamples( size_t num_trIN, size_t num_baIN )    
   {    
      std::vector<double>* ptr_v   ;
      for (size_t t1 = 0 ; t1 < num_trIN ; t1++) // for each treatment
      {
        std::list<std::vector<double>* > T_l_of_v  ;    // create the list to store the vectors (_set_ba.size())
        for (size_t t2 = 0 ; t2 < num_baIN ; t2++) // create the batch vectors
        {
           ptr_v = new std::vector<double> ;  // dynamic allocation of a vector<double>
           T_l_of_v.push_back(ptr_v) ;        // push the pointer to the empty vector into the list 
        }
          
        this->_llpvd_T.push_back(T_l_of_v) ; // push the list of vectors onto the main list
      }
   }


	 /*! return the pointer to the vector of doubles corresponding to the wanted_treatment from the wanted_batch */
	 std::vector<double> * Get_VectorsOfSamples( size_t wanted_treatmentIN, size_t wanted_batchIN ) 
	 {
			std::list<std::vector<double>* > T_l_of_v ;
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists     
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;

	    it_T=_llpvd_T.begin() ;
			std::advance(it_T,wanted_treatmentIN) ;
			T_l_of_v = *it_T ;
			it_T_l_of_v=T_l_of_v.begin() ;
			std::advance(it_T_l_of_v,wanted_batchIN) ;
			return (*it_T_l_of_v) ;
	 }

	/*!
   * \details Deletes/Frees the vectors<double> arrays that hold scores data. Used by destructor.
   * \see Return_T_WithVariableData()
   */
  void Free_T_VectorsOfVariableData( void )
  {
    // look through the _llpvd_T list of lists and deallocate the vector<double> data in them.
      std::vector<double>* t_v   ;
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists     
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      for (it_T=this->_llpvd_T.begin(); it_T!=this->_llpvd_T.end(); ++it_T) // for each treatment
      {
        std::list<std::vector<double>* > T_l_of_v = *it_T  ;    // 
        for (it_T_l_of_v=T_l_of_v.begin(); it_T_l_of_v!=T_l_of_v.end(); ++it_T_l_of_v) // iterate through the batches vectors
        {
            t_v  = *it_T_l_of_v ;
            delete t_v ;
        }
      }
  }

 /*!
   * \details Zeros the lengths of the vectors<double> arrays that hold scores data, so we can fill them with another set of Scores data
   * \see Return_T_WithVariableData()
   */
  void Empty_T_VectorsOfScores( void )
  {
    // look through the _llpvd_T list of lists and deallocate the vector<double> data in them.
      std::vector<double>* v   ;
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists     
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      for (it_T=this->_llpvd_T.begin(); it_T!=this->_llpvd_T.end(); ++it_T) // for each treatment
      {
        std::list<std::vector<double>* > T_l_of_v = *it_T  ;    // 
        for (it_T_l_of_v=T_l_of_v.begin(); it_T_l_of_v!=T_l_of_v.end(); ++it_T_l_of_v) // iterate through the batches vectors
        {
            v  = *it_T_l_of_v ;
            v->clear() ;
        }
      }
  }
} ;  // end of class CExperimentData



  
/*!
   * Creates data structures to store treatment and batch information with helper information _TB_array[][] and TB_treatment_info[]
	 * and mapping functions to get non-0 indexed treatment and batch numbers into a 0 based format. Also stores total combinations of
	 * treatments. To be used in conjunction with a CExperimentData object, that holds a the data (variable or PC) that is used in the simulateion
	 * of distribution functions.
   *  \brief Creates data structures to store treatment and batch information.
   */
class  CExperimentStructure {
private:
  
  size_t *  _TB_array_block ;    /*! Stores the 1D block of data that _TB_array[][] will use to create a 2D version*/
  size_t ** _TB_array ;          /*! _TB_array[][] is a 2D array of numbers of samples present in the _llpvd_T structure 
                                      rows are treatments, cols are batches, and values at [tr][ba]
                                      are number of times this treatment 'tr' was present in batch 'ba' */  
  
  double *  _noofcombs_block ;    /*! the bulk memory to store the 2D array data for _noofcombs[][] 
                                      \see Create_noofcombs_Array() */
  double ** _noofcombs ;           /*! binomial theorm nchoosek() results - The number of combinations possible for each set of treatment(s) in a batch selected out 
                                      of set of all common treatments. Values are arranged [tr][ba] = treatments in rows and batches in columns (this is a transposed version of what is used in Matlab code) 
                                      \see Create_noofcombs_Array() 
																			\verbatim 
																			Example	for a treatment/batch structure = 4 treatments and 3 batches
																								_TB_array[tr][ba] =				
																																	  ba1, b12, ba3 | _TB_treatment_count[tr]  (== row sums)
																															tr1 :		{{2,	2,	2},		|			6
																															tr2 :		 {2,	2,	2},		|			6
																															tr3 :		 {4,	2,	2},		|			8
																															tr4 :		 {0,	2,	2}}		|			4
																												             --------------
																								_TB_batch_count[ba] =	 {8,	8,	8}
																								(== column sums)

																								_noofcombs[tr][ba] =	(N,K) = binomial coefficient of set of size N, choose K	
																								where N == _TB_treatment_count[tr] and K == TB[tr][ba]
																																		   ba1,			b12,		 ba3
																															tr1 :	{{(6,2),		(6,2),	(6,2)},
																															tr2 :	 {(6,2),		(6,2),	(6,2)},
																															tr3 :	 {(8,4),		(8,2),	(8,2)},
																															tr4 :	 {(4,0),		(4,2),	(4,2)	}}				
																			\endverbatim
																			*/  


  double *  _totalcombs ;							/*! Total number of combinations of each available set of treatments in a batch (== column products of _noofcombs[][tr] 
																					\see Create_noofcombs_Array() 
																			*/  
  
  size_t * _TB_treatment_count ;				/*! TB_treatment_info[] is the row sums of _TB_array[][] - each value of which
																						equals the total number of a particular treatment in all batches */
  size_t * _TB_batch_count ;						/*! _TB_batch_count[] is the column sums of _TB_array[][] - each value of which
																						equals the total number samples in a particular batch  
																			*/
  

 protected:  
	 // Copy of input data
  std::vector<double> v_treatmentinfo ;		/*! Holds information about a samples treatment	- one value for each expected sample within the experiment */
  std::vector<double> v_batchinfo ;		/*! Holds information about a samples batch			- one value for each expected sample within the experiment */
  size_t _num_samples ;						/*! The total number of PCA samples of the study input data */

	// derived information
  std::set<size_t> _set_tr ;		/*! Set of unique treatments */
  std::set<size_t> _set_ba ;		/*! Set of unique batchs */
	size_t _num_tr ;							/*! The number of unique treatments, equivalent to _set_tr.size() */
  size_t _num_ba ;							/*! The number of unique batches, equivalent to _set_ba.size() */
  
  // map of the unique elements to range 0..(set.size()-1)
  std::map<size_t, size_t> _map_tr ; /*! Map of unique treatments to sequential values 0,1,2,..   */
  std::map<size_t, size_t> _map_ba ; /*! Map of unique batchs to sequential values 0,1,2,..  */
  
  

  
  /*!
   * \details Releases the assigned memory from the _TB_array[][]
   * \see Return_T_WithVariableData()
   */
  void Free_TRBA_Array( void ) 
  {
     if (this->_TB_array != NULL)
        delete [] _TB_array ;
     if (this->_TB_array_block != NULL)
        delete [] _TB_array_block ;
     if (this->_TB_treatment_count != NULL)
        delete [] _TB_treatment_count ;
     if (this->_TB_batch_count != NULL)
        delete [] _TB_batch_count ;
      if (this->_noofcombs_block  != NULL)
        delete [] _noofcombs_block  ;
      if (this->_noofcombs != NULL)
        delete [] _noofcombs ;
      if (this->_totalcombs != NULL)
        delete [] _totalcombs ;
  }
  
 
  

  /*! Creates the (empty) arrays that are later populated with batch/treatment size values.
   * * _noofcombs[][] is a 2D array of numbers of combinations of choosing K treatments from a set of N (common) treatments present in the study.
     * Currently, the nchoosek() result is calculated in CSimulateBatchDistribution() using the CFactorial class \see class CSimulateBatchDistribution
     * _totalcombs[] is a vector of products of columns of _noofcombs[][], which is the total number of possible combinations of all treatments selected
     * 'tr' ast a time in a particular batch, where 'tr' is the number of treatments of the particular type of interest (==K).
\verbatim     
Create_noofcombs_Array() creates:  
				_noofcombs[][] = {{NA, NA, NA},  
												 {NA, NA, NA},  
												 {NA, NA, NA},
												 {NA, NA, NA}}

				_totalcombs[] = { NA, NA, NA }

And when filled with experment structure data:
				TR [][] = {{2, 2, 2},
                   {2, 2, 2},
                   {4, 2, 2},
                   {0, 2, 2}}
      
      then _noofcombs[][] = {{15, 15, 15},  
                            {15, 15, 15},  
                            {70, 28, 28},
                            {0,  6,  6 }}
       (N.B binmial coef. (N,K) : (6,2) = 15; (8,4) = 70 ; (8,2) = 28 ; (4,2) = 6) )
      
      and _totalcombs[] = { 15 x 15 x 70, 15 x 15 x 28 x 6, 15 x 15 x 28 x 6 } = { 15750, 37800, 37800 }


\endverbatim
   */
    void Create_noofcombs_Array( void )
    {
        // all these arrays are created here, as called from the class creator method destroyed in the destructor and filled with data in class CSimulateBatchDistribution
        this->_noofcombs_block =  (double *) new double[this->_num_tr * this->_num_ba] ;  // allocate block that will fit whole 2D array
        this->_noofcombs = (double **) new double*[this->_num_tr] ;									// assign a vector of pointers - one for each row (treatment type)
        for(size_t i = 0; i < this->_num_tr; i++)                            // loop over the treatments
          this->_noofcombs[i] = (double*) this->_noofcombs_block  + (this->_num_ba * i) ;  // set pointers to offsets into bulk memory area
         
       this->_totalcombs  =  (double *) new double[this->_num_ba] ;						// create the total combs array
    }

  
  
public:
  
  // Accessor functions
  size_t            GetNumTreatments( void )  {   return  this->_num_tr ; }       /*! Accessor function - returns the total number number of treatments in the experiment */
  size_t            GetNumBatches( void )     {   return  this->_num_ba ;  }      /*! Accessor function - returns the total number number of batches in the experiment */
  size_t            GetNumSamples( void )     {   return  this->_num_samples ; }  /*! Accessor function - returns the total number of samples that are present in the experiment */
  std::vector<double> * GetTreatmentInfoVect( void ) {   return &this->v_treatmentinfo ; }	/*! Accessor function - returns the vector of information about a samples treatment	- one value for each expected sample within the experiment */
  std::vector<double> * GetBatchInfoVect( void )    {   return &this->v_batchinfo ; }		/*! Accessor function - returns the vector of information about a samples batch			- one value for each expected sample within the experiment */
  std::map<size_t, size_t>* GetTreatmentMap( void ) {   return &this->_map_tr ;}  /*! Map of unique treatments to sequential values 0,1,2,..   */
  std::map<size_t, size_t>* GetBatchMap( void ) {   return &this->_map_ba ;}  /*! Map of unique batchs to sequential values 0,1,2,..  */
  
  /*! Accessor function - returns how many samples are in the batch of interst /param batchIN 0 based batch of interest */
	size_t            GetBatcheSizes(size_t batchIN)			  
	{  
    if (this->_TB_batch_count != NULL) 
    return  this->_TB_batch_count[batchIN] ; 
		else return 0 ;
  }
  /*! Accessor function - returns how many samples are of a particular treament /param treatmentIN 0 based treatment type of interest */
	size_t            GetTreatmentSizes(size_t treatmentIN)   
	{ 
    if (this->_TB_treatment_count != NULL)  
    return  this->_TB_treatment_count[treatmentIN] ; 
		else return 0 ;
  }
	size_t **         GetTB_array( void )        {   return this->_TB_array ;  }		 /*! 2D array containing sample numbers for a particular treatment (along rows) and batch (down columns) */
  double **         GetNoofcombsArray( void )  {   return this->_noofcombs ;  }   /*! Accessor function - returns 2D array containing the number of combinations possible for each set of treatment(s) in a batch selected out of the set of all common treatments  */
	double *          GetTotalcombsArray( void ) {   return this->_totalcombs ;  }  /*! Accessor functions - Total number of combinations of each available set of treatments in a batch (== column products of _noofcombs[][tr] */
 
 
   /*!
   * \details Default constructor.
   */
  CExperimentStructure()
  {
      _TB_array = NULL ;
      _TB_treatment_count = NULL ;
      _TB_batch_count = NULL ;
  }
  
  /*! **Constructor method**
   * \details // constructor that provides all required data -
   *    \param std::vector<size_t> tr : Vector of integer numbers describing a samples treatment value. Each sample has a treatment and batch number as they belong to a batch and haveundergone a specific treatment.
	 *		\param std::vector<size_t> ba : Vector of integer numbers describing a samples batch number. Each sample has a treatment and batch number.
   *    This information is used to create a T object which is a list of list of vectors, with each vector containing the PCA scores for 
   *    a specific treatment in a specific batch (of one PC of interest at a time). If a batch did not have an occurence of a particular treatment then the vector exists, 
   *    however it contains no scores data. (\see Return_T_WithVariableData for an example of the structure of the "list of lists of pointers to vectors of doubles" 
	 *		which could be regarded as a 'jagged' 2D array).
   *    Selection of all samples with a common treatment is done by selecting (from _llpvd_T) the inner <list of vectors> from the outer list item corresponding
   *    to the treatment of interest. This selection is akin to selecting a 'row' of a 2D matrix.  Mapping of treatment numbers to the outer list row is done with the 
   *    std::map<> functionality. This is required as treatment and batch numbers do not necesarily start at 0, as required for access of C/C++ array objects.
   *    \param std::vector<double> vect_scores : Is a full set of scores, of dimension [ number of samples x number of PCs] and its layout specified by dataOriented
   *    \param dataOriented. The vect_score data is copied to a CPCAScoresArray (\see class PCAScoresArray) object that copies the data and has access methods 
   *    that allow the selection of specific PCs.
   * \result On creation the obect fills the _llpvd_T structure with data from PC1 of vect_scores<>
   
	\verbatim    
	 Data
	 {SamplesAsRows=0, SamplesAsCols=1}
    SamplesAsRows: (S = sample, v = variable (==PC))
          V1  V2  V3  V4
      S1  -   -   -   -
      S2  -   -   -   - 
      S3  -   -   -   -
  \endverbatim    
   */
  CExperimentStructure(std::vector<size_t> &trIN, std::vector<size_t>& baIN, size_t samplesIN)
  {
			assert(trIN.size() == baIN.size()) ;

      this->_TB_array = NULL ;
      this->_TB_treatment_count = NULL ;
      this->_TB_batch_count = NULL ;
      this->_TB_array_block = NULL ;
      this->_noofcombs_block  = NULL ;
      this->_noofcombs  = NULL ;
      this->_totalcombs = NULL ;
      

      // Copy input data
      this->v_treatmentinfo.assign(&trIN[0], &trIN[0]+trIN.size());
      this->v_batchinfo.assign(&baIN[0], &baIN[0]+baIN.size());
      
      this->_num_samples = samplesIN ;		/*!  We will assert that any input vector containing data will be of this std::vector::size() */
      this->DetermineUniqueTrBa(trIN, baIN ) ;	// this will give us the number of treatments and number of batches
      this->MapInputValuesToStandardRange() ;		// re-maps non-0 based and non-sequentialy numbered treatment and batch numbers back to start at 0
      	     
			this->Create_noofcombs_Array() ;					// creates the array but does not populate with values
  }
  
  
  ~CExperimentStructure()
  {
      Free_TRBA_Array() ;
      this->_TB_array = NULL ;
      this->_TB_array_block = NULL ;
      this->_TB_treatment_count = NULL ;
      this->_TB_batch_count = NULL ;
      this->_noofcombs_block  = NULL ;
      this->_noofcombs  = NULL ;
      this->_totalcombs = NULL ;
  }
  

 
 /*!
   * Determines how many unique treatments and batches there are from the input vectors
   */
  void DetermineUniqueTrBa(std::vector<size_t>& v_trIN, std::vector<size_t>& v_baIN )
  {
    for (size_t i=0; i< v_trIN.size(); i++) 
        _set_tr.insert(v_trIN[i]);    // creates a set of unique values only (if insertion value is a duplicate then it is ignored)
        
    for (size_t i=0; i< v_baIN.size(); i++) 
        _set_ba.insert(v_baIN[i]);    // creates a set of unique values only (if insertion value is a duplicate then it is ignored)
   
    this->_num_tr = _set_tr.size() ; // _set_tr.size() are the number of unique treatments
    this->_num_ba = _set_ba.size() ; // _set_ba.size() are the number of unique batches
  }
  
  /*!
   * \details Maps possibly disjoint and variable range treatment and batch information to continous range that starts at 0 (seq_val)
   */
  void MapInputValuesToStandardRange ( void )
  {
      assert(v_treatmentinfo.size() == v_batchinfo.size()) ;  /*! \assert Each sample must have batch and treatment information. So, the vectors of treatments and batches should be equal in size */
			std::set<size_t>::iterator it1 ;

      size_t seq_val = 0 ;  // set this as 0 as we use the mappings to determin C/C++ array offsets which start at 0.
      for (it1=_set_tr.begin(); it1!=_set_tr.end(); ++(it1))  // we know that *it values are unique as they are from a <set>
      {
         size_t tint = *it1 ;
         this->_map_tr.insert(std::pair<size_t,size_t>(tint, seq_val));    // map unique values to a sequence starting at 1
         seq_val++ ;
      }
      seq_val = 0 ;
      for (it1=_set_ba.begin(); it1!=_set_ba.end(); ++(it1)) // we know that *it values are unique as they are from a <set>
      {
         size_t tint = *it1 ;
          this->_map_ba.insert(std::pair<size_t,size_t>(tint, seq_val));    // map unique values to a sequence starting at 1
         seq_val++ ;
      }  
  }
  
   /*!
   * \details Creates:
   * _TB_array[][] is a 2D array of the numbers of values present in the _llpvd_T vector<double> structures (i.e. the sizes of vectors of doubles)
   *      rows are treatments, cols are batches, and values at [tr][ba]
   *      are number of times this treatment 'tr' was present in batch 'ba'
   * _TB_treatment_count[] 1D array of integers which store the counts of values in _TB_array[][] rows
   *      so is the number of samples having a particular treatment that are present in all batches.
   * _TB_batch_count[] 1D array of batch sizes
	 * 
	 * \requires	Requires the data to be added to the _llpvd_T structure using Return_T_WithVariableData(), otherwise 
	 *						it fills the TB array with 0's
   */
  void Create_TB_Array(CExperimentData * llpvd_T_IN )
  {
      this->_TB_array_block =  (size_t *) new size_t[_num_tr * _num_ba] ;  // allocate block that will fit whole 2D array
      
      _TB_array = (size_t **)new size_t*[_num_tr] ;  // assign a vector of pointers - one for each row (treatment type)
      for(size_t i = 0; i < _num_tr; i++)  // loop over the treatments
        _TB_array[i] = (size_t *) _TB_array_block  + (_num_ba * i) ;  // set pointers to offsets into bulk memory area
        
        
      std::list<std::vector<double>* > T_l_of_v ;
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists     
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      // _TB_array stores the lengths of the vectors in the list of list of vectors (_llpvd_T)
      size_t tr_i = 0 ; 
      size_t ba_i = 0 ;
      for(it_T=llpvd_T_IN->ReturnData()->begin(); it_T != llpvd_T_IN->ReturnData()->end(); ++it_T)  // iterate over the (outer) treatment list - to return a std::list
      {
         T_l_of_v = *it_T ;     // this is a std::list<std::vector<double>* >
         ba_i = 0 ;
         for(it_T_l_of_v=T_l_of_v.begin(); it_T_l_of_v != T_l_of_v.end(); ++it_T_l_of_v)   // iterate over each inner (common batch) list - to return a std::vector<double> of sample data with common treatment and batch
         {
            size_t v_size = (*it_T_l_of_v)->size() ; // the size of the vector 
            _TB_array[tr_i][ba_i] = v_size ;       // place size in the 2D array
            ba_i++ ;  // increment batch index
          }
        tr_i++ ;  // increment treatment index
      }
      
      _TB_treatment_count  = (size_t*) new size_t[_num_tr ] ;
      _TB_batch_count      = (size_t*) new size_t[_num_ba ] ;
			// initialise arrays to zero
			for(tr_i = 0; tr_i < _num_tr; tr_i++) this->_TB_treatment_count[tr_i] =  0 ;
			for(ba_i = 0; ba_i < _num_ba; ba_i++) this->_TB_batch_count[ba_i] = 0 ; 
      
			for(tr_i = 0; tr_i < _num_tr; tr_i++)  // loop over the treatments
      {
        for(ba_i = 0; ba_i < _num_ba; ba_i++)  // loop over batches
        {
            this->_TB_treatment_count[tr_i] += _TB_array[tr_i][ba_i] ; // number of a particular treatment type that are present in all batches
            this->_TB_batch_count[ba_i] += _TB_array[tr_i][ba_i] ;
        }  
      }
            
  }  // end of method Create_TB_Array() 

   
    /*!
   * \details Returns all samples that have an identical treatment
   * \note "Return Value Optimisation" should optimise the return structure without having to duplicate it.
      R test code:    
          TestHarmanBTExtraction(0.95, 1, 0, -1) # balanced data and the first treatment. 
          # Output:  [1]  -0.5303412 -14.3886991  -2.3057925  -7.6598890 -61.2716444 -40.7700399 
          TestHarmanBTExtraction(0.95, 2, 1, -1) # unbalanced data and the second treatment.
          # Output: [1]   0.336662  29.410803   7.054404   8.726244 -16.731484 -10.862511
          TestHarmanBTExtraction(0.95, 2, 2, -1) # anotherunbalanced data set and the third treatment. 
          # Output: [1]   9.463589   7.264324  56.980556  22.915228 -16.469382 -12.538409 -12.994833   2.660779
          TestHarmanBTExtraction(0.95, 2, 3, -1) # anotherunbalanced data set and the fourth treatment. 
          # Output: [1]  26.094083 28.784813 -1.564443 -1.604018
        
   */
   std::vector<double> GetCommonTreatmentsSamples(CExperimentData* llpvd_T_IN , size_t treatmentIN, bool mappingNeeded  )    
   {
      std::vector<double> ret_vect ;
      
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists
      std::list<std::vector<double>* > T_l_of_v ;
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      
      size_t tr_index ; 
      if (mappingNeeded == true)
        tr_index = _map_tr[treatmentIN] ;
      else
        tr_index = treatmentIN ;
        
      assert(tr_index < this->_num_tr ) ;  // tr_index is 0 based
      
			it_T = llpvd_T_IN->ReturnData()->begin();         // assign an iterator to the start of the outer list.
      std::advance (it_T,tr_index);         // iterate to the treatment of interest
      T_l_of_v = *it_T ;                    // assign the list 
      it_T_l_of_v = T_l_of_v.begin();       // obtain an iterator for the (inner) list variable
      ret_vect.reserve(ret_vect.size() + _TB_treatment_count[tr_index]);  // use _TB_treatment_count[] to pre-allocate all the space needed
        
      for (size_t t1 = 0 ; t1 < this->_num_ba ; t1++) // for each batch
      {          
          std::vector<double>* vec = *it_T_l_of_v ; // obtain the vector that relates to this sample          
          //ret_vect.reserve(ret_vect.size() + distance(vec->begin(),vec->end()));  // reserve extra space in the result vector (could use _TB_treatment_count[] to pre-allocate all the space needed)
          ret_vect.insert(ret_vect.end(),vec->begin(),vec->end());        
          std::advance (it_T_l_of_v,  1);    // advance the iterator to move to the next vector relating to this treatments batch
      }
          
     return ret_vect ;

   }
   
 
	 
	 /*! overloaded version that will extract the batch data of a single batch from a given _llpvd_T structure and place it in a single std::vector<double> 
	    Data is placed in vector in ascending treatment order */
	 std::vector<double> GetCommonBatchSamples(CExperimentData* llpvd_T_IN , size_t batchIN, bool mappingNeeded  )    
   {
			std::vector<double> ret_vect ;
      
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists
      std::list<std::vector<double>* > T_l_of_v ;
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      
      size_t ba_index ;
      if (mappingNeeded == true)
        ba_index = this->_map_ba[batchIN] ;
      else
        ba_index = batchIN ;
        
      assert(ba_index < this->_num_ba ) ;  // ba_index is 0 based
   
			it_T = llpvd_T_IN->ReturnData()->begin();                       // assign an iterator to the start of the outer list.
			for (size_t t1 = 0 ; t1 < this->_num_tr ; t1++)  // for each treatment
      {  
          T_l_of_v = *it_T ;
          it_T_l_of_v = T_l_of_v.begin();           // obtain an iterator for the (inner) list variable
          std::advance (it_T_l_of_v,ba_index);      // use the iterator to move to the vector relating to this treatments batch
          std::vector<double>* vec = *it_T_l_of_v ; // obtain the vector that relates to this batch and treatment
          
          ret_vect.reserve(ret_vect.size() + distance(vec->begin(),vec->end()));
          ret_vect.insert(ret_vect.end(),vec->begin(),vec->end());
          std::advance (it_T,1);         // iterate to the next list of treatments
      }
          
     return ret_vect ;

   }

 
  /*! An object factory to create CExperimentData objects from input vectors of scores data
   * \details Populates the CExperimentData objects '_llpvd_T' structure - a list of list of pointers to vectors of doubles
	 *				The structure is described as follows:
	 *					- The outer list is a list of <lists of vectors<double>*>, 
	 *					- There is one one <lists of vectors<double>*> for each treatment type. 
	 *					- each <lists of vectors<double>*> contains a vector<double>*, one for each batch
	 *					- each  vector<double>* is populated by the PCA 'scores' data for samples that corresponds to the treatement type and batch in which it was processed
	 *
   *    e.g.  So, for 4 treatments occuring in 3 batches: 
   *					e.g.  t = treatment b = batch
              list[ list( t1_b1<2.1,3.2>,						t1_b2<3.2,0.8>,			t1_b3<0.1,3.1> ),      
                    list( t2_b1<1.8,1.8>,						t2_b2<2.1,-3.5>,		t2_b3<-2.2,1.2> ),   
                    list( t3_b1<2.5,3.4,2.1,-3.5>,	t3_b2<5.1,-2.2>,		t3_b3<3.2,-0.8> ),
                    list( t4_b1<>,									t4_b2<2.1,3.2>,			t4_b3<3.2,-2.8> ) ]
   */
	CExperimentData* Return_T_WithVariableData(std::vector<double>* vect_scoresIN)
  {
    // look through the tr and ba vectors, translate to the mapped integer sequence
      size_t tr_key, ba_key ;
      size_t tr_index, ba_index ; 
      
			CExperimentData * llpvd_T_OUT = (CExperimentData *) new  CExperimentData( this->GetNumTreatments(), this->GetNumBatches() ) ;

      assert(this->_num_samples == vect_scoresIN->size()) ; /*! \assert The input vector of scores must be the same size as the number of samples expected
                                                            (i.e. one PCA score for each sample) */      
      std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists
      std::list<std::vector<double>* > T_l_of_v ;
      std::list<std::vector<double>* >::iterator it_T_l_of_v ;
      for (size_t v1 = 0; v1 < this->_num_samples; v1++) // for each sample
      {
        tr_key = this->v_treatmentinfo[v1] ;        // get the map key from the input dataset (which is sample number v1 treatment)
        ba_key = this->v_batchinfo[v1] ;        // get the map key from the input dataset (which is sample number v1 batch)
        tr_index = _map_tr[tr_key] ;         // the mapped index into the treatment list
        ba_index = _map_ba[ba_key] ;         // the mapped index into the batch list
				it_T = llpvd_T_OUT->ReturnData()->begin();				// assign an iterator to the start of the list of lists
        std::advance (it_T,tr_index);       // use the iterator to move to the list that relates to this treatment
        T_l_of_v = *it_T ;                  // assign the list to a referenceable list variable
        it_T_l_of_v = T_l_of_v.begin();     // obtain an iterator for the list variable
        std::advance (it_T_l_of_v,ba_index);// use the iterator to move to the vector relating to this treatments batch
        std::vector<double>* vec = *it_T_l_of_v ; // obtain the vector that relates to this sample
        vec->push_back((*vect_scoresIN)[v1]) ;  // place the PCA score data into the vector
		//		_PRINTSTD << "tr_key: " << tr_key << " tr_index " << tr_index << " ba_key: " << ba_key << " ba_index " << ba_index <<  " value: " << vect_scoresIN->data()[v1]  <<  std::endl ;
      }
       //_PRINTSTD << std::flush ;

			return llpvd_T_OUT ;
  }
	 


} ;  // end of class CExperimentStructure definition



/*! This class collects together two objects - 
			- A pointer to an CExperimentStructure object (that can be reused in multiple CExperimentWithPCAData objects)
      - And a CExperimentData object that holds the data structure to hold a single variable (in this case a single PC) of an experiments multivariate data.
				This object holds a list of list of pointers to vectors of doubles (llpvd)  structure.
				For implementation, this class requires the help of a CPCAScoresArray object to feed in the variable data 
		It allows reuse of the Expermint structure information with different data so that we can process multiple objects  in parallel using multiple CPU threads. 
		\brief Combines a CExperimentStructure and CExperimentData to define the structure of the experiment and provide a single PC's worth of data for simulating the data distribution.
		*/
class  CExperimentWithPCAData 
{
private :     
	CExperimentStructure	* ExptStruct ;
	CExperimentData				* ExptData ;		/*! _llpvd_T structure with data */
	size_t								_PCnum ;				/*! keep this incase we want to reference what PC is being processed */
public:   
	CExperimentStructure	* GetExptStructureObject() { return this->ExptStruct ; } 
	CExperimentData	* GetExptDataObject()  { return this->ExptData ;}  
	size_t GetPCNum()  { return this->_PCnum ;}  
	size_t            GetNumTreatments()  {   return  this->ExptStruct->GetNumTreatments()  ; }       /*! Accessor function - returns the total number number of treatments in the experiment */
	size_t            GetNumBatches()     {   return  this->ExptStruct->GetNumBatches() ;  }      /*! Accessor function - returns the total number number of batches in the experiment */
	size_t            GetNumSamples()     {   return  this->ExptStruct->GetNumSamples() ; }  /*! Accessor function - returns the total number of samples that are present in the experiment */
  

	/* fill the _llpvd_T structure with data from input vector */
	CExperimentData	* SetExptDataObjectData(std::vector<double> * pcascoresIN)  
	{ 
		delete this->ExptData ;
		this->ExptData   = this->ExptStruct->Return_T_WithVariableData(pcascoresIN) ;
    
    return ( this->ExptData ) ;
	} 


		/*! **Constructor method** */
		CExperimentWithPCAData(CExperimentStructure	* ExptStructIN, std::vector<double> * pcascoresIN, size_t fillWithPCNumIN ) 			
		{
			this->_PCnum = fillWithPCNumIN ;
			if (ExptStructIN != NULL)
				this->ExptStruct = ExptStructIN ;
			else
			{
				_PRINTERROR << "CExperimentWithPCAData::constructor error: CExperimentStructure cannot be NULL" << std::endl ;
			//	error( ) ; // exit(-1) ;
			}
			// add some scores data to this object
			this->ExptData   = this->ExptStruct->Return_T_WithVariableData(pcascoresIN) ;  // This fills the _llpvd_T structure with data from input data
			ExptStruct->Create_TB_Array(this->ExptData) ;  // Created last as TB requires a correct _llpvd_T structure that describes the numbers of samples and treatments in each batch.      
		}

		/*! Destructor Method */
  ~CExperimentWithPCAData()
  {
		delete this->ExptData ;
  }


	
	/*! 
	* SetExptDataFromBatchVectors() This method places the input data back into the CExperimentData llpvd structure 
	* The input is a pointer to a vector of vectors of doubles that is created using the CMatrixData::CreateCMatrixFromCExperimentData() 
	* which uses CExperimentStructure::GetCommonBatchSamples() method	
	* \see CMatrixData::CreateCMatrixFromCExperimentData()
	* \see CExperimentStructure::GetCommonBatchSamples()
	*/ 
	void SetExptDataFromBatchVectors( std::vector<std::vector<double> > * vvd_matrix_of_batchscoresIN ) 
	{
			assert(vvd_matrix_of_batchscoresIN->size() == ExptStruct->GetNumBatches() ) ;  // there should be one inner vector of doubles for each batch present in the experiment

			std::vector<std::vector<double> >::iterator vd_it_batchscores ;
			vd_it_batchscores = vvd_matrix_of_batchscoresIN->begin() ;
			
			size_t curr_element ;
			for (size_t i_ba = 0 ; i_ba < ExptStruct->GetNumBatches() ; i_ba++)
			{
							
				std::vector<double> * ptr_vd = &*vd_it_batchscores  ;
				double * ptr_input_data = ptr_vd->data() ; // access to the array of doubles from the batch
				curr_element = 0 ;
				for (size_t i_tr = 0 ; i_tr < ExptStruct->GetNumTreatments() ; i_tr++)
				{
						std::vector<double>* vec = this->ExptData->Get_VectorsOfSamples(i_tr,i_ba) ; // obtain the vector that relates to this batch and treatment
						size_t numTreatments ;
						numTreatments = vec->size() ;

						double * ptr_output_data = vec->data() ;
						for (size_t t_cpy =0 ; t_cpy < numTreatments; t_cpy++ ) 
							ptr_output_data[t_cpy] = ptr_input_data[curr_element++] ;
				}
				std::advance (vd_it_batchscores,  1); 				
			}
	}
	
	/*! 
	 *  Method to copy over CExperimentWithPCAData scores values corresponding to unique batches into a vvd vector
	 *  This method has an inverse SetExptDataFromBatchVectors() that populates a CExperimentWithPCAData->CExperimentData llpvd array with data from a vvd array
	 */
	std::vector< std::vector<double> > * CreatePVVD_FromCExperimentData( void )
  {
		std::vector<std::vector<double> > * ret_pvvd = (std::vector<std::vector<double> > *) new std::vector<std::vector<double> >() ;

		for (size_t batch_int =0 ; batch_int < this->GetNumBatches(); batch_int++)
    {
			ret_pvvd->push_back( this->GetExptStructureObject()->GetCommonBatchSamples(this->GetExptDataObject(), batch_int, false)) ;
    }  

		return ret_pvvd ;
/*
    for (size_t batch_int =0 ; batch_int < ptr_ExptIN->GetNumBatches() ; batch_int++)
    {
      CalculateMeanAndStore( this->vvd_matrix[batch_int], batch_int  ) ; 
    }*/
  }
  
  
  /*! Maps back from the vectors of samples with common batch and treatment values to the original input vector of data sample order
   * \details Uses the CExperimentData objects '_llpvd_T' structure - a list of list of pointers to vectors, the input order vectors
   * of batch and treatment and the maping structures to get the input vectors rferenced back to 0 
   *    e.g. _llpvd_T for 4 treatments occuring in 3 batches: 
   *	 t = treatment b = batch
   list[ list( t1_b1<2.1,3.2>,						t1_b2<3.2,0.8>,			t1_b3<0.1,3.1> ),      
         list( t2_b1<1.8,1.8>,						t2_b2<2.1,-3.5>,		t2_b3<-2.2,1.2> ),   
  list( t3_b1<2.5,3.4,2.1,-3.5>,	t3_b2<5.1,-2.2>,		t3_b3<3.2,-0.8> ),
  list( t4_b1<>,									t4_b2<2.1,3.2>,			t4_b3<3.2,-2.8> ) ]
  */
  std::vector<double> * Return_vector_from_T( void  )
  {
    // look through the tr and ba vectors, translate to the mapped integer sequence
    size_t tr_key, ba_key ;
    size_t tr_index, ba_index ; 
    size_t ba_tr_count ;
     
    std::vector<size_t> * v_ba_tr_Counts = (std::vector<size_t> *) new  std::vector<size_t>( this->ExptStruct->GetNumTreatments()  * this->ExptStruct->GetNumBatches(), 0 ) ;
    std::vector<double> * pvd_samples_OUT = (std::vector<double> *) new  std::vector<double>( this->ExptStruct->GetNumSamples() ) ;
    std::list< std::list<std::vector<double>* > > * llpvd_T_ptr = this->ExptData->ReturnData() ;  // CExperimentData
    
    std::list< std::list<std::vector<double>* > >::iterator it_T ;  // Assign an iterator to the start of the list of lists
    std::list<std::vector<double>* > T_l_of_v ;
    std::list<std::vector<double>* >::iterator it_T_l_of_v ;
    for (size_t v1 = 0; v1 <  this->ExptStruct->GetNumSamples(); v1++) // for each sample
    {
      tr_key = this->ExptStruct->GetTreatmentInfoVect()->at(v1) ;        // get the map key from the input dataset (which is sample number v1 treatment)
      ba_key = this->ExptStruct->GetBatchInfoVect()->at(v1)  ;        // get the map key from the input dataset (which is sample number v1 batch)
      tr_index = this->ExptStruct->GetTreatmentMap()->find(tr_key)->second ;         // the mapped index into the treatment list - find() returns an iterator
      ba_index = this->ExptStruct->GetBatchMap()->find(ba_key)->second ;         // the mapped index into the batch list - find() returns an iterator
      
      ba_tr_count = v_ba_tr_Counts->at(ba_index * this->ExptStruct->GetNumTreatments() + tr_index ) ;
      v_ba_tr_Counts->at(ba_index * this->ExptStruct->GetNumTreatments() + tr_index )++ ; 
      
      it_T = llpvd_T_ptr->begin();				// assign an iterator to the start of the list of lists
      std::advance (it_T,tr_index);       // use the iterator to move to the list that relates to this treatment
      T_l_of_v = *it_T ;                  // assign the list to a referenceable list variable
      it_T_l_of_v = T_l_of_v.begin();     // obtain an iterator for the list variable
      std::advance (it_T_l_of_v,ba_index);// use the iterator to move to the vector relating to this treatments batch
      std::vector<double>* vec = *it_T_l_of_v ; // obtain the vector that relates to this sample
      pvd_samples_OUT->at(v1) = ((*vec)[ba_tr_count]) ;  // place the score data into the return vector
      //		_PRINTSTD << "tr_key: " << tr_key << " tr_index " << tr_index << " ba_key: " << ba_key << " ba_index " << ba_index <<  " value: " << vect_scoresIN->data()[v1]  <<  std::endl ;
    }
  //_PRINTSTD << std::flush ;
  delete v_ba_tr_Counts ;
    
  return pvd_samples_OUT ;
  }

} ;


// N.B. These to classes (CFactorial and CMapSelectKFromN) are include from file CMapSelectKFromN.h
/* \brief Class that computes the factorial of positive input integers and stores them as an array for later access
 *  \reference 
 * 
 */
//template <class T> class CFactorial // From include file CMapSelectKFromN.h

// typedef std::pair<size_t, size_t> key_NK_pair ;  // requires #include <utility>
  
/* \brief Creates a map of all recursively defined binomial sequences relating to choosing K values from the set of N.
 * 
 *  \details Creates all binomial sequences for input numbers (N,K), N being the number in a set and K being the size of the sample to take form the set.
 *    Each sequence is the scan of the binomial coefficients (N-1,K-1),(N-2,K-1),(N-3,K-1)...(N-X,K-1) where X = N-K-1
 *    e.g. (8,4)  =  (7,3),(6,3),(5,3),(4,3),(3,3) where X = 8-4-1 = 3
 *                =   35  , 20,   10,   4,    1
 *          scan  =   35,   55,   65,  69,  70
 *    Results in a map of each sequence scan from N..1 and K..1
 * 
 *  \export
 */
// class CMapSelectKFromN // From include file CMapSelectKFromN.h



/*! \breif Equivalent of a std::vector<double>* accompanied by the number of values used to make up a single double in the std::vector<> so we can calculate means correctly later.
 * Used to hold the (equivalent of the) a single batch of the X matrix - one CVectorSum for each treatment present in a batch, along with _numSums, the number of values that the sum is made up of.
* _numSums value is needed so we know how many values are added so we can calculate a correct mean in the M array
* \param _numSums How many values were summed to make each value in the vector
*/
class CVectorSum 
{

private: 
	std::vector<double>* v_sums ;  /*! each element is the sum of _numSums values, which equals the number of treatments of a particular kind present in a batch */
	size_t _numSums ;   /*! this is how many values were summed to make each value in the vector */

public: 
	CVectorSum(size_t numSumsIN)
	{
		this->_numSums = numSumsIN ;
		this->v_sums = (std::vector<double> *) new std::vector<double> ;   // dynamic allocation of a vector<double>
	}
	~CVectorSum() 
	{
		delete this->v_sums ;
	}
	/*! accessor method */
	size_t GetNumSums( void ) { return _numSums ; } 
	/*! accessor method */
	std::vector<double>* GetV_SumPtr( void ) { return this->v_sums ; } 
	
	/*! assign method */
	void  SetNumSums( size_t numSumsIN ) { this->_numSums = numSumsIN ; } 
	
  /*! assign method */
	 void AssignV_SumPtr( std::vector<double>* newVect ) 
	 {  
		 delete this->v_sums ;  // remove the old data
		 this->v_sums = newVect ; // point to new input data
	 } 

} ;



/*! 
 \brief Calculates the average, variance and standard deviation of input values
 * Class to compute running means and variances
 * Modified from a website that took it from Knuth TAOCP vol 2, 3rd edition, page 232
 */
class CRunningStat
{
	private:
        size_t _count ;  /*! the number of values currently input to the calculation */
        double _oldM, _newM ; /*! means */
				double _oldS, _newS;  /*! variances */

   public:
        CRunningStat() : _count(0) {}

    void Clear()
    {
        _count = 0;
    }

      void Push(double input)
      {
          _count++;

          // See Knuth TAOCP vol 2, 3rd edition, page 232
          if (_count == 1)
          {
              _oldM = _newM = input;
              _oldS = 0.0;
          }
          else
          {
              _newM = _oldM + (input - _oldM)/_count;
              _newS = _oldS + (input - _oldM)*(input - _newM);
    
              // set up for next iteration
              _oldM = _newM; 
              _oldS = _newS;
          }
      }

      size_t NumDataValues() const
      {
          return _count;
      }

      double Mean() const
      {
          return (_count > 0) ? _newM : 0.0;
      }

      double Variance() const
      {
          return ( (_count > 1) ? _newS/(_count - 1) : 0.0 );
      }

      double StandardDeviation() const
      {
          return sqrt( Variance() );
      } 

			void operator=(const CRunningStat& other) //  assignment operator
			{		
					/* copy data from other's storage to this storage */
					this->_count	= other._count				;
					this->_oldM		= other._oldM		;
					this->_newM		= other._newM		;
					this->_oldS		= other._oldS		;
					this->_newS		= other._newS		;										
			}
};  // end of class CRunningStat



/*! \brief Class that computes the averages of the combinatorial range of a set of numbers - Simulates the distribution of an input set of numbers to return a simulated average and standard deviation.
 *   This class is currently used by the CProcessScores class (which also creates the CExperimentStructure object that is passed into the creator function )
 *   MATLAB file simbat_unbalancedYO.m
 * 
 */
class  CSimulateBatchDistribution
{
private:

	// This object does not need to be destroyed (in this class)
	CExperimentWithPCAData * ptr_ExptData ;						/*! Pointer to a pre-created CExperimentWithPCAData object - used for the TB_ptr[][] information, nocombs[][] and _totalcombs[] */

	// All following objects need to be created and destroyed
	CFactorials<double> * _factn ;											/*! CFactorial object will be created in constructor.  */
	size_t _arrayAlign ;															/*! Array length multiple - so array boundaries are set to good values for vectorisation (used in CMapSelectKFromN::GetScan_Vectorised() ) */
	std::vector<CMapSelectKFromN *> vect_NKMap ;			/*! The CMapSelectKFromN object lets us access array indecies of a binomial sequence */
	
	// These are used for calculation and storage of means and standard deviations of the M data
	std::vector<CRunningStat *> v_M_mean_and_stddev ;	/*! This is the end goal of this class - the mean and standard deviations of the column of M, one for each batch. */
	std::vector<double> v_batchmeans ;     /*! One value for each batch - filled using GetVectsOfMeansAndStddevs() */
	std::vector<double> v_batchstdevs ;    /*! One value for each batch - filled using GetVectsOfMeansAndStddevs() */

	// Original Matlab code named data structures
	std::vector< CVectorSum * > v_X ;  /*! This is the  X matrix in original Matlab code, although the entries are sums of treatment scores, not the each individual treatment score (that would be summed later when creating the M matrix results) */
	std::vector<double> v_M ;		   /*! This is the (row averages of the) M matrix in original Matlab code. Stores all combinations of the v_X vector data */

	CSelectRandom<double> * _RandomSelect ;

  
	bool bool_FillCombsCalled ; /*  */

	/*! Calculates the mean and standard deviation of the the M[] array data
	*  \note Only to be called from CreateMMatrix() 
	*/
	void CalculateRunningStats( void )
	{
			double * M_ptr = v_M.data() ;
			CRunningStat * ptr_RS_stats ;
			ptr_RS_stats = new CRunningStat() ;
			for (size_t seq_val = 0 ; seq_val < this->v_M.size() ; seq_val++ )  // for each value in the loop
			{
				ptr_RS_stats->Push(M_ptr[seq_val]) ;
			}
			this->v_M_mean_and_stddev.push_back(ptr_RS_stats) ;
	}

  /*! Clear the data from the objects v_M_mean_and_stddev vector of CRunningStat * objects
	*/
	void FreeRunningStatsObjects( void )
	{
			for (size_t i_rs = 0 ; i_rs < this->v_M_mean_and_stddev.size() ; i_rs++ )  // for each value in the loop
			{
				if (v_M_mean_and_stddev[i_rs] != NULL) 
					delete v_M_mean_and_stddev[i_rs] ;
			}

			this->v_M_mean_and_stddev.clear() ;
	}


 /*!  Check if a N,K map already exists in the vector of current CMapSelectKFromN maps 
 *  \param nIN the N value to check for in the CMapSelectKFromN * object
 *	\param kIN the K value to check for in the CMapSelectKFromN * object
 */
	CMapSelectKFromN * GetNKMap(size_t nIN, size_t kIN)
	{
		CMapSelectKFromN * t_nkMap ;
		CMapSelectKFromN * t_nkMap_ret ;
		t_nkMap_ret = NULL ;
	//	MakeN_LTE_K
		 for (std::vector<CMapSelectKFromN *>::iterator it=this->vect_NKMap.begin(); it!=this->vect_NKMap.end(); ++it)
		 {
			   t_nkMap =  *it ;  // look at each CMapSelectKFromN in the vector
				 bool n_bool = (t_nkMap->GetN() == nIN) ;
				 bool k_bool = (t_nkMap->GetK() == kIN) ;

				 if ((n_bool == true) && (k_bool == true))
					t_nkMap_ret = t_nkMap ;  // return the matching CMapSelectKFromN pointer
		 }
		 return (t_nkMap_ret) ;
	}

	
	/*!
	*  Copies a previously calculated CRunningStats object that was placed in the vector v_M_mean_and_stddev
	*  and places the copy onto v_M_mean_and_stddev vector
	*  \info Used when a short circuit for a calulation is determined when a previous treatment structure is identical to a current one
	*   so we can copy the results of the previous one
	*/
	CRunningStat * CopyRunningStats( size_t batchToCopy )
	{
			CRunningStat * ptr_RS_stats ;
			CRunningStat * ptr_RS_stats_copy ;

			ptr_RS_stats_copy = new CRunningStat() ;  // create a new CRunningStat object
			ptr_RS_stats = this->v_M_mean_and_stddev[batchToCopy]  ;  // obtain the CRunningStat object to copy
			*ptr_RS_stats_copy = *ptr_RS_stats ;  // Copy contents of most recent into the new copy version
			this->v_M_mean_and_stddev.push_back(ptr_RS_stats_copy) ;  // place new copy into the vector<>
			
			return (ptr_RS_stats) ;
	}

	CRunningStat * GetRunningStats( size_t batch )
	{
			CRunningStat * ptr_RS_stats ;			
			ptr_RS_stats = this->v_M_mean_and_stddev[batch]  ;
			return (ptr_RS_stats) ;
	}

	double ReturnRunningStatsAverage( size_t batch )
	{
			CRunningStat * ptr_RS_stats ;			
			ptr_RS_stats = this->v_M_mean_and_stddev[batch]  ;
			return (ptr_RS_stats->Mean()) ;
	}
	double ReturnRunningStatsStdDevs( size_t batch )
	{
			CRunningStat * ptr_RS_stats ;			
			ptr_RS_stats = this->v_M_mean_and_stddev[batch]  ;
			return (ptr_RS_stats->StandardDeviation()) ;
	}

	/*! Each CVectorSum store how many values were added to form each item in its vector<double> object. 
	 * 	This function adds these values from each CVectorSum present to be used as the 'row average' divisor
	 */
	double GetCurrentDivisor( void ) 
	{
			CVectorSum * v_s ;
			double current_divisor ;
			current_divisor = 0.0 ;
			for (std::vector<CVectorSum *>::iterator it=this->v_X.begin(); it!=this->v_X.end(); ++it)
  		{
					v_s = *it ;
					current_divisor += v_s->GetNumSums() ;
			}
			return current_divisor ;
	}

	/*! checks that the input double is within the range of values that can be stored by type size_t */
	size_t GetTotalAsSize_t( double d_total_cocombs ) 
	{
		  size_t t_total_cocombs ;

			if (d_total_cocombs >  ((size_t) -1) )
			{
#if (_PRINT==0)
				 _PRINTSTD << "Error: CSimulateBatchDistribution::GetTotalAsSize_t(): Number of permutations exceeds the capacity of size_t data type" << std::endl;
#endif
				// TODO: This state implies we need to invoke a different random selection proceedure than the one currently provided.
				// exit(-1) ;
				t_total_cocombs =	((size_t) -1) ;
			}
			else
				t_total_cocombs = (size_t) d_total_cocombs ;

			return t_total_cocombs ;
	}
 

	/*!   Set up the X vector information - Used in creation of the M matrix
	 *		- Create the number of combinations on addition of each X CVectorSum of data
	 *		- Place the vector of double into a 2D array structure for easy access
	 */ 
	void FillXVectSizeAndPtrs(size_t *  XVectSizes_ptr, double ** XVectors_ptrs	)
	{
		CVectorSum * v_s ;
		size_t i ;
		i = 0 ;
		for (std::vector<CVectorSum *>::iterator it=this->v_X.begin(); it!=this->v_X.end(); ++it)
		{
  			v_s = *it ;
				// place the pointers to the first elements of the vector<double> X matrix data in the array of double *
				// XVectors_ptrs[i] = &(v_s->GetV_SumPtr()->begin()[0]) ; // that is pretty ugly!
				XVectors_ptrs[i]		= v_s->GetV_SumPtr()->data() ;
				if (i == 0)
					XVectSizes_ptr[i] = v_s->GetV_SumPtr()->size() ;  // this is the progression so we can calculate the array index from a single index 
				else // if (i > 0)
					XVectSizes_ptr[i] = v_s->GetV_SumPtr()->size() * XVectSizes_ptr[i-1] ;
				i++ ;
		}
	}

	/*! Checks if the batch structure of a previous batch is identical
	\return The batch which has identical structure or the current batch if it is so far unique */
	size_t CheckIfColumnIndentical(size_t ba_num)
	{
		size_t ** TB_ptr  = this->ptr_ExptData->GetExptStructureObject()->GetTB_array() ;  // this holds all the treatments and batch numbers (each column is a batch and each row is a common treatment)      
			size_t numTr			= this->ptr_ExptData->GetExptStructureObject()->GetNumTreatments() ;
		//	size_t maxBatch		= this->ptr_ExptData->GetExptStructureObject()->GetNumBatches() ;
			bool		isSameAsPrevious  ;
			size_t	batchToCopy ;

			batchToCopy = ba_num ;  // initialise the result to equal the input batch.
		//	assert(ba_num < maxBatch) ; // make sure the batch number of interest is within the range of the experiment structure data
			
			// These loops check if column values are all identical for all columns before the batch column of interest
			for (size_t ba = 0 ; ba < ba_num ; ba++)  // for each column that is before the (batch) column of interest
			{
				isSameAsPrevious = true  ;
				for (size_t tr = 0 ; tr < numTr ; tr++) // for each row in the column
				{
							isSameAsPrevious &= (TB_ptr[tr][ba] == TB_ptr[tr][ba_num]) ; // check if the values they contain are the same and update boolean result
				}
				if (isSameAsPrevious == true)
					batchToCopy = ba ;
			}

			return (batchToCopy) ; // N.B. if the only column that matches is the batch of interest itself then we have to calculate the X and M data
	}


public: 

	/*! Creator - accepts the experiment batch and treatment description object
	 * 
	 * \param arrayAlignIN - does not do much unless we use the 'vectorised' binomial calculation function
	 * \param pointer to the experiment batch and treatment structure
	 * \param seedIN - if 0, the uses the current time as a seed - low precision clock used, so all instances may get the same seed.
	 */
  CSimulateBatchDistribution(size_t arrayAlignIN,  CExperimentWithPCAData * ptr_Expt_IN, size_t seedIN ) 
  {
    this->_arrayAlign = arrayAlignIN ;
    this->_factn = (CFactorials<double> *) new CFactorials<double>(20) ;
    ptr_ExptData = ptr_Expt_IN  ;
		this->bool_FillCombsCalled = false ;

	this->_RandomSelect = (CSelectRandom<double> * ) new CSelectRandom<double>(seedIN) ;

  }

	/*! destructor
	*/
	~CSimulateBatchDistribution() 
	{
		delete this->_factn ;
		CVectorSum * temp_vs ;
		// delete the batch vectors - the X matrix may be empty as it is cleared at the end of method CreateMMatrix()
		for (std::vector<CVectorSum *>::iterator it=this->v_X.begin(); it!=this->v_X.end(); ++it)
		{
			temp_vs = *it ;
			delete temp_vs ; // this will call destructor of CVectorSum objects so will delete the CVectorSum->v_sums<> vector data
		}

		this->v_X.clear() ;
		this->v_M.clear() ;

		CRunningStat * ptr_RS_stats ;
		for (std::vector<CRunningStat *>::iterator it=this->v_M_mean_and_stddev.begin(); it!=this->v_M_mean_and_stddev.end(); ++it)
		{
			  ptr_RS_stats = *it ;
				delete ptr_RS_stats ; // this will call destructor of CVectorSum objects so will delete the CVectorSum->v_sums<> vector data
		}
		this->v_M_mean_and_stddev.clear() ; // clear the array RunningStats objects used to calculate means and stddev's

		CMapSelectKFromN * t_NKmap ;
		for (std::vector<CMapSelectKFromN *>::iterator it=this->vect_NKMap.begin(); it!=this->vect_NKMap.end(); ++it)
		{
			  t_NKmap =  *it ;
				delete t_NKmap ;				
		}
		this->vect_NKMap.clear() ;

		delete this->_RandomSelect ;

	}  // end of destructor


	CExperimentWithPCAData * GetCExptPCAPointer( void ) 
	{
		return (this->ptr_ExptData) ;
	}

	/*! This does most of the work - simulates each batch to obtain the simulated batch mean and stddev*/
	void SimulateBatches(size_t max_combinations, bool forceRandomIN) 
	{

		this->FreeRunningStatsObjects() ;
			// for each batch
		for ( size_t ba_num = 0 ; ba_num < ptr_ExptData->GetExptStructureObject()->GetNumBatches() ; ba_num++ )
		{
			this->CalculateOrCopySimulatedMeans(ba_num, max_combinations, forceRandomIN) ;
			//this->CreateXMatrix(ba_num, max_combinations) ;
			//this->CreateMMatrix(ba_num, max_combinations) ;
		}
		this->ComputeVectsOfMeansAndStddevs() ;
		//this->PrintMeans() ;
		//this->PrintStdDevs() ;
	}

	 /*! \details Calaculates the number of combinations possible (e.g. the binomial N choose K) 
	 *					for the number of treatments present in a batch (K) from the total number of that treatment type (N)
	 *					and stores these values in the 2D array nocombs_ptr[][], along with the total combinations (for each batch) in array totalnocombs_ptr[][]
	 *					Also creates the vector of CMapSelectKFromN * pointers to all N,K sequence maps
	 *  \param nIN the N value to check for in the CMapSelectKFromN * object
	 *	\param kIN the K value to check for in the CMapSelectKFromN * object
	 */
	 void FillCombinationsMatrix( size_t max_combinations )
  {
      size_t ** TB_ptr            = ptr_ExptData->GetExptStructureObject()->GetTB_array() ;  // this holds all the treatments and batch numbers
      double ** nocombs_ptr       = ptr_ExptData->GetExptStructureObject()->GetNoofcombsArray() ; // we will fill a column ([][ba]) of this array each time this function is called
      double *  totalnocombs_ptr  = ptr_ExptData->GetExptStructureObject()->GetTotalcombsArray() ;// we will form a single element of this array ([ba]) each time this function is called
      double    d_nocombs ;
      //size_t    t_nocombs ;
     
			CMapSelectKFromN * temp_newNKMap ;
			
			for (size_t ba=0; ba < ptr_ExptData->GetExptStructureObject()->GetNumBatches() ; ba++) 
			{
				totalnocombs_ptr[ba]  = 1.0 ;  // initilise this to 1.0 as we multiply each value iteratively
				for (size_t tr=0; tr < ptr_ExptData->GetExptStructureObject()->GetNumTreatments() ; tr++) // create the batch vectors
				{
					size_t n = ptr_ExptData->GetExptStructureObject()->GetTreatmentSizes(tr) ;   // total number of treatments of type t1 
					// _TB_treatment_count[] 1D array of integers which store the counts of values in _TB_array[][] rows
					size_t k = TB_ptr[tr][ba]  ;            // treatments of type t1 in this batch

					temp_newNKMap = this->GetNKMap(n,k) ; // check if the (N,K) pair already exists in the vector of maps
					if (temp_newNKMap == NULL) // if the N,K pair does not exist, then create a new CMapSelectKFromN object and place it in the vector for reuse if possible
					{
						// The CMapSelectKFromN object lets us access array indices of a binomial sequence
						// We may want to only create this once and add N,K pairs to the map as different ones are needed ( N,K == 300,7 took about 0.3 seconds to create)
						CMapSelectKFromN * temp_newNKMap = (CMapSelectKFromN *) new CMapSelectKFromN(std::make_pair(n,k), this->_arrayAlign) ;  // this object lets us access array indecies of a binomial sequence
						vect_NKMap.push_back(temp_newNKMap) ;
					}

					d_nocombs  = this->_factn->nchoosek(n,k) ;  // this could possibly be very large - if so invoke a random selection routine. Or could be 1
   
					assert(d_nocombs <=  (size_t) -1) ; // '(size_t) -1' ends up the largest unsigned integer of size sizeof(size_t)
					if (d_nocombs >  (size_t) -1) 
					{
#if (_PRINT==3)					  
						_PRINTERROR << "Error: CMapSelectKFromN::FillCombinationsMatrix(): Number of combinations exceeds the capacity of size_t data type" << std::endl;
#endif
					}
				
					nocombs_ptr[tr][ba]   = d_nocombs ;
					totalnocombs_ptr[ba]  *= d_nocombs ;
				}
			}

			this->bool_FillCombsCalled = true ;
  }
  
  /*! \brief Creates (the equivalent of) the X matrix (from original Matlab code), which is the stepping point to create the M matrix of all combinations
   *  \detals Stores a vector of doubles for each treatment type that is present in the batch of interest, and each vector 
	 *  contains a set of sums computed from values taken in lots of the number of treatments present (k) in the batch of 
	 *  interest from the set of all possible samples of that treatment type (n)
   *  \reference MATLAB file simbat_unbalancedYO.m
   */
  void CreateXMatrix( size_t batch_num, size_t max_combinations, bool forceRandomIN )
  {
      //assert(this->bool_FillCombsCalled == true) ;
	  if (this->bool_FillCombsCalled == false)
		this->FillCombinationsMatrix(max_combinations) ;  /// creates arrays nocombs[][] and total_combs[][] and the vector of CMapSelectKFromN * pointers to all N,K sequence maps
		
	  CVectorSum * ptr_v_s  ;    // a pointer to a std::vector<double> * and a int stating how many values each vect[] item actually is the sum of.

	  size_t ** TB_ptr            = ptr_ExptData->GetExptStructureObject()->GetTB_array() ;  // this holds all the treatments and batch numbers
    double ** nocombs_ptr       = ptr_ExptData->GetExptStructureObject()->GetNoofcombsArray() ; // we will fill a column ([][ba]) of this array each time this function is called
    size_t  t_nocombs ;
    double   d_test_nocombs ;
    CMapSelectKFromN * temp_NKMap ; // The CMapSelectKFromN object lets us access array indecies of a binomial sequence
		bool b_MAPPINGNEEDED ;				// We explicitly use the mapped range of 0...NumTreatments-1 so no need to map from original (possibly disjoint) treatment numbers.
		b_MAPPINGNEEDED = false ;
		for (size_t tr=0; tr < ptr_ExptData->GetExptStructureObject()->GetNumTreatments() ; tr++) // create the batch vectors
		{
			size_t n = ptr_ExptData->GetExptStructureObject()->GetTreatmentSizes(tr) ;   // n is the total number of treatments of type tr
			// _TB_treatment_count[] 1D array of integers which store the counts of values in _TB_array[][] rows
			size_t k = TB_ptr[tr][batch_num]  ;      // treatments of type tr in this batch
			ptr_v_s = (CVectorSum *) new CVectorSum(k) ;  // this stores how many values were summed to make each value in the vector

			// We may want to only create this once and add N,K pairs to the map as different ones are needed ( N,K == 300,7 took about 0.3 seconds to create)
			//  newNKMap = (CMapSelectKFromN *) new CMapSelectKFromN(std::make_pair(n,k), this->_arrayAlign) ;  // this object lets us access array indecies of a binomial sequence
			temp_NKMap = this->GetNKMap(n,k) ;
			double current_sum ;

			d_test_nocombs = nocombs_ptr[tr][batch_num]   ;  // these values have passed an assert() that they are within the size possible for a size_t value in the FillCombinationsMatrix() function
			t_nocombs = this->GetTotalAsSize_t(d_test_nocombs)  ;
			std::vector<double> * v_vectd ;
			//std::vector<size_t> * v_random64bitRangeIN  ;
#if (_PRINT==1)
			_PRINTSTD << "CreateXMatrix() t_nocombs: " << t_nocombs << " \t max_combinations: " << max_combinations << std::endl  ;
#endif          
           
		if ((t_nocombs > max_combinations) | (forceRandomIN==true)) 
		{
#if (_PRINT==1)
		  _PRINTSTD << "CreateXMatrix() entered forceRand=T section " << std::endl  ;
#endif
          
            // Get a vector of common treatment samples - returns a vector of sample values from common treatments (treatment == tr) from each batch
						
            std::vector<double> v_treatments = ptr_ExptData->GetExptStructureObject()->GetCommonTreatmentsSamples( ptr_ExptData->GetExptDataObject(), tr, b_MAPPINGNEEDED  ) ; 
#if (_PRINT==1)
          _PRINTSTD << "CreateXMatrix() obtained v_treatments vector " << std::endl  ;
#endif
            // Get the random indexes for use to choose the treatment data and create averages 
            // This is not very efficient - using 64 bits to store an index that possibly will not be much greater than 500
            // Some optimisation has been done
            if (k != 0)
            {
				      v_vectd = _RandomSelect->SelectWITHOUTREPLACEMENTReturnRowSumFast(&v_treatments,  k,  max_combinations ) ;
            }
            else
            {
            	v_vectd  = (std::vector<double> *) new std::vector<double>(0) ;
            }
			      ptr_v_s->AssignV_SumPtr(v_vectd) ;  // if k is zero (none of a particular treatment is present in this particular batch) then this will be an empty vector.
#if (_PRINT==1)
            _PRINTSTD << "CreateXMatrix() finished forceRand=T section " << std::endl  ;
#endif
		 
		 }
     else  // do full combinatorial version
     {
    			ptr_v_s->GetV_SumPtr()->reserve(t_nocombs) ;
    			ptr_v_s->GetV_SumPtr()->assign(t_nocombs,0.0) ;
    			double * ptr_Xdata  = ptr_v_s->GetV_SumPtr()->data() ;
    			// Get a vector of common treatment samples
    			std::vector<double> v_treatments = ptr_ExptData->GetExptStructureObject()->GetCommonTreatmentsSamples( ptr_ExptData->GetExptDataObject(), tr, b_MAPPINGNEEDED  ) ; // returns a vector of sample values from common treatments (treatment == tr) from each batch
    			std::vector< size_t > v_index(k) ; // one for each treatment to be choosen (size() == k)

#if (_PRINT==3)   			
    			double tr_sum  ;
    			tr_sum  = 0.0 ;
					_PRINTSTD << "v_M average for ba = " << batch_num << " tr = " << tr << " : " ; 
    			for (size_t t1 = 0 ; t1 < v_treatments.size() ; t1++)
					{
    			  tr_sum += v_treatments[t1] ;
						_PRINTSTD << v_treatments[t1] << ", "  ;
					}
    			if (v_treatments.size() > 0)
    			  _PRINTSTD << "  ave: " <<  (tr_sum / (double) v_treatments.size()) << std::endl  ;
#endif   			
    			
    			for (size_t seq_val = 1 ; seq_val <= t_nocombs ; seq_val++)
    			{
      				// requires a sequence that starts at 1
      				// This returns a vector of array indexes
      			//	v_index = temp_NKMap->GetScan_shortcut(seq_val) ;  // this copies the index data so we can run the loop in parallel
      			v_index = temp_NKMap->GetScan_shortcut_debug(seq_val, (std::string)" CreateXMatrix() ") ;
      				current_sum = 0.0 ;
      				for (size_t sum = 0 ; sum < v_index.size() ; sum++)
      				{
      					  current_sum += v_treatments.at(v_index.at(sum)) ;
      				}
      				//ptr_v_s->GetV_SumPtr()->push_back(current_sum) ;
							v_index.clear() ;
      				ptr_Xdata[seq_val-1] = current_sum ;
    			}
		
         }   // end of full combinatorial version
         v_X.push_back(ptr_v_s) ;  // push the vectorSum into the vector/array. One for each treatment (could be 
      }
      
  }
  

  /*! \brief Creates the row averaged version of the M matrix, by sampling from the X matrix, either completly or via random selection if number of combinations is too high.
   *  \details Stores a vector of doubles for each treatment type that is present in the batch of interest, and each vector contains a set of averages computed from
   *  values taken in lots of the number of treatments present (k) in the batch of interest from the set of all possible 
   *  samples of that treatment type (n)
   *  \reference MATLAB file simbat_unbalancedYO.m
   *  \note Only needs to be called if the treatment structure is different to all previously calculated batches
   */
  void CreateMMatrix( size_t batch_num, size_t max_combinations, bool forceRandomIN )
  {
		CVectorSum * v_s  ;    		// a temporary pointer so we can iterate through a vector o them I understand a little better about what is happening
		double ** nocombs_ptr       = ptr_ExptData->GetExptStructureObject()->GetNoofcombsArray() ; // a column ([][ba]) for each batch
		double *  totalnocombs_ptr  = ptr_ExptData->GetExptStructureObject()->GetTotalcombsArray() ;// one for each batch

		double d_total_cocombs ;
		size_t t_total_cocombs ;

		// allocate working values and get pointers to the data arrays of the X matrix
		size_t *  XVectSizes_ptr  = (size_t *) new size_t[v_X.size()] ;		// this is the progression so we can calculate the array index from a single index
		double ** XVectors_ptrs		= (double **) new double*[v_X.size()]  ; // assign a vector of pointers to the X data vectors of doubles - one for each row (treatment type)
		double *  M_ptr ;  // this will store the resulting M data for this batch
		// size_t *  XVectSizes_ptr, double ** XVectors_ptrs

		// DO NOT USE FillXVectSizeAndPtrs() FOR RANDOM SELECTION - the Sizes are multiples of previous and current for > [0][]
		this->FillXVectSizeAndPtrs(XVectSizes_ptr, XVectors_ptrs ) ;  // fills the 2D array with X array data (X vector data can have different lengths so this is an equivalent of a 'jagged' array)

		//std::vector<double> * v_vectd ;
		//std::vector<size_t> * t_v_random64bit ;
		std::vector<unsigned int> *  t_v_random32bit ;

		// these values have passed an assert() that they are within the size possible for a size_t value in the FillCombinationsMatrix() function
		d_total_cocombs = totalnocombs_ptr[batch_num]   ;
		t_total_cocombs = GetTotalAsSize_t(d_total_cocombs) ;

		double current_divisor ; // used to create the mean of the data - precalculate it here as each M value is made up of the same number of component X values
		//size_t numXVects = 0 ;
		current_divisor = this->GetCurrentDivisor() ;  // sums all the constituent X CVectorSum numSums values

		
		if ((t_total_cocombs > max_combinations) | (forceRandomIN==true)) 
		{

				//	_PRINTSTD << "CSimulateBatchDistribution::CreateMMatrix() Message: t_total_cocombs is large, so we will invoke a random selection of values to simulate the distribution. t_total_cocombs= " <<  t_total_cocombs<< std::endl ;
				this->v_M.reserve(max_combinations) ;
				this->v_M.assign(max_combinations,0.0) ;

				M_ptr = this->v_M.data() ;
				CVectorSum * temp_vs ;
				size_t i_XvectSize ;

				// This loop is parallelisable
				// for each value in the loop, create indexes into the jagged X array (XVectors_ptrs[][]),
				// and average the data
				// for the remaining X vectors, we use integer division to access the correct element
				for (size_t t1 = 0 ; t1 < this->v_X.size()  ; t1++)  // for each X vector  - I was using numXVects. 1 vector for each treatment
				{
					// I should test if it is more efficient to create the array of indexes or create them individually and discard them (with no memory bandwidth hit)
					// could also test for required index size (e.g. byte, short, int of 64 bit int) and create appropriate indexes of that size.
					// Simple test of size_t vs unsigned int shows *no difference* in speed
					temp_vs = v_X[t1] ;
					i_XvectSize = temp_vs->GetV_SumPtr()->size() ;
					// Get the random indexes for use to choose the treatment data and create averages

					//t_v_random64bit = _RandomSelect->ReturnVectOf64bitIntegersInRange_CSTDLIB(max_combinations, 0, 0, i_XvectSize  ) ;
					t_v_random32bit = _RandomSelect->ReturnVectOf32bitIntegersInRange_CSTDLIB(max_combinations, 0, 0, i_XvectSize  ) ;
					// this will return NULL if the 0, i_XvectSize range is 0 => no data to select from(or -ve)
					//if (t_v_random64bit != NULL)
					if (t_v_random32bit != NULL)  // if NULL then i_XvectSize is probably = 0
					{
						// size_t * ptr_Xindex =  t_v_random64bit->data() ;

						unsigned int * ptr_Xindex =  t_v_random32bit->data() ;
						for (unsigned int  seq_val = 0 ; seq_val < max_combinations ; seq_val++ )
						{
								M_ptr[seq_val] += XVectors_ptrs[t1][ptr_Xindex[seq_val]] ; // yes, using an array variable full of random indexes is a performance sin
						}
						// delete t_v_random64bit ;
						delete t_v_random32bit ;
					}
				}

				// now create the equivalent 'row mean' of the M matrix from the Matlab code
				for (size_t seq_val = 0 ; seq_val < max_combinations ; seq_val++ )
				{
					M_ptr[seq_val] /= current_divisor ;
				}
				
			}
			  else  // do full combinatorial version - select
			  {
				this->v_M.reserve(t_total_cocombs) ;
				this->v_M.assign(t_total_cocombs,0.0) ;
				M_ptr = this->v_M.data() ;
				// This loop is paralelisable
				// for each value in the loop, create indexes into the jagged X array (XVectors_ptrs[][]),
				// and average the data
				for (size_t seq_val = 0 ; seq_val < t_total_cocombs ; seq_val++ )   
				{
					double current_sum ; // used to store the temporary M data results
					size_t t_Xindex ;    // used to index into the X vectors

					// the first vector index (==[0][]) is the 'remainder' so use modulo division
					t_Xindex = seq_val % XVectSizes_ptr[0] ;
					current_sum =	XVectors_ptrs[0][t_Xindex] ;
					// for the remaining X vectors, we use integer division to access the correct element
					for (size_t t1 = 1 ; t1 < this->v_X.size()  ; t1++)  // for each X vector  - I was using numXVects
					{
						// t_Xindex = seq_val / XVectSizes_ptr[t1] ;  // integer division should result in the correct rounding (floor())
						t_Xindex = (seq_val  / XVectSizes_ptr[t1-1]) % ((size_t) nocombs_ptr[t1][batch_num]) ;
						current_sum += XVectors_ptrs[t1][t_Xindex] ;
					}
					// now create the equivalent 'row mean' of the M matrix from the Matlab code
					M_ptr[seq_val] = current_sum / current_divisor ;
				}
		}   // end of full combinatorial version
  

	// Now calculate the average and stddev of the M vector data.
	this->CalculateRunningStats() ;

		// tidy up
    // delete the X matrix data vectors
     for (std::vector<CVectorSum *>::iterator it=this->v_X.begin(); it!=this->v_X.end(); ++it)
  	 {
			  v_s = *it ;
			  delete v_s ;
		 }
		 this->v_X.clear() ;
		 this->v_M.clear() ; // this will be reused for each batch
		 delete [] XVectSizes_ptr ;
		 delete [] XVectors_ptrs  ; // but do not delete the memory that the pointers point to
  }  // end of CreateMMatrix()
  


	/*!
	 *  This method either calculates the simulated mean and stddev or copies it from a previously computed batch if it is an identically structured batch.
	 */
	void CalculateOrCopySimulatedMeans(size_t ba_num, size_t max_combinations, bool forceRandomIN )
	{

			size_t batchToCopy ;				
			bool weHaveToCalaculate ;
			// N.B. if the only column that matches is the batch of interest itself (==ba_num) then we have to calculate the X and M data
			batchToCopy = this->CheckIfColumnIndentical(ba_num) ; // returns ba_num if no other columns match 

			weHaveToCalaculate = (batchToCopy == ba_num) ;

			if (weHaveToCalaculate == true)
			{
				this->CreateXMatrix(ba_num, max_combinations, forceRandomIN) ;
				this->CreateMMatrix(ba_num, max_combinations, forceRandomIN) ;
			}
			else
			{
				this->CopyRunningStats(batchToCopy) ;
			}
	}

	/*! 
	 *  This fills the member variable vectors with the means and stdevs of the data found in v_M. 
	 *  This is equivalent to the average of the rows of the Matlab 'M' array, being averaged and standard deviation calcuated.
	 */
	void ComputeVectsOfMeansAndStddevs( void )
	{
		this->v_batchmeans.clear() ;
		this->v_batchstdevs.clear() ;

		for (size_t t1 = 0 ; t1 < this->v_M_mean_and_stddev.size() ; t1++)
		{
			this->v_batchmeans.push_back(ReturnRunningStatsAverage(t1)) ;
			this->v_batchstdevs.push_back(ReturnRunningStatsStdDevs(t1)) ;
		}
	}

  void PrintMeans( void ) 
	{
		_PRINTSTD << "v_M column means:  " ;
		for (size_t t1 = 0 ; t1 < this->v_batchmeans.size(); t1++)
			_PRINTSTD << this->v_batchmeans[t1] << " " ;
		_PRINTSTD << std::endl ;
	}

	void PrintStdDevs( void ) 
	{
		_PRINTSTD << "v_M column stddevs:  " ;
		for (size_t t1 = 0 ; t1 < this->v_batchmeans.size(); t1++)
			_PRINTSTD << this->v_batchstdevs[t1] << " " ;
		_PRINTSTD << std::endl ;
	}

	std::vector<double> * GetPtrMeans( void ) 
	{
		return &this->v_batchmeans ;
	}

	std::vector<double> * GetPtrStdDevs( void ) 
	{
		 return &this->v_batchstdevs ;
	}
} ; // end class CSimulateBatchDistribution





/*! \brief Helper class to calculate means of batches - stores all common batch samples in a single vector so we can calculate the batch means quickly.
  Interoperates with the CExperimentWithPCAData object to reorganise the _llpvd_T data into the batch only data.
 */
class  CMatrixWithMeans
{
  public: 
  std::vector<std::vector<double> > * pvvd_matrix ; /*! Each inner vector are the scores of all samples in a batch (from any treatment) */
  std::vector<double> vd_means ; /*! The average of the inner std::vectors of vvd_matrix */
  double _scale ; /*! The scale factor used on the means of the vvd_matrix vectors as used in the optimistaion stage of Harman */
  bool _meanCentred ; /*! If true then the current means vector has been removed from the data */
  
  /*! **Constructor method**
   */
  CMatrixWithMeans( void )
  {
		this->pvvd_matrix = NULL ;
    this->_scale = 0.0 ;
    this->_meanCentred = false ;
  }
  
  /*! **Constructor method**
   * 
   * \param std::vector<double> input_scores - returned from GetCommonBatchSamples( size_t batch, bool mappingNeeded )  
	 * \param CExperimentStructure * ptr_Expt - contains the structure of the experiment
   * \warning Requires input_scores that have been ordered via a CExperimentStructure object. Do not use with the 'raw' input scores ( use CExperimentStructure::GetCommonBatchSamples() ).
   * \see CExperimentStructure::GetCommonBatchSamples() as the method for ordering the scores
  */
  CMatrixWithMeans(std::vector<double>& input_scores, CExperimentStructure * ptr_Expt)
  {
		this->pvvd_matrix = NULL ;
    this->_scale = 0.0 ;
    this->_meanCentred = false ;

    this->pvvd_matrix = (std::vector<std::vector<double> > *) new std::vector<std::vector<double> >() ; 
   
		if (input_scores.size() == ptr_Expt->GetNumSamples())
    {
    
      // ret_vect.reserve(ret_vect.size() + distance(vec->begin(),vec->end()));
      std::vector<double>::iterator input_scores_it ;
      input_scores_it = input_scores.begin() ;
                
      for (size_t i =0 ; i < ptr_Expt->GetNumBatches() ; i++)  // for each batch present, create a vector containing the samples it holds
      {
        size_t batchSize = ptr_Expt->GetBatcheSizes(i) ;
        std::vector<double> temp(batchSize) ; 
        temp.insert(temp.begin(),input_scores_it,input_scores_it+batchSize); // copy section of std::vector<double> input data
        this->pvvd_matrix->push_back(temp) ;
        // calculate the sum and the mean
        this->CalculateMean(temp,i) ;
        //double s = std::accumulate(temp.begin(),temp.end(), 0.0);
        //this->vd_means[i] = s / (double) temp.size() ;
        std::advance (input_scores_it,batchSize ) ;
      }
    }
  }

	/*! **Destructor method**
   */
  ~CMatrixWithMeans( void )
  {
    this->_scale = 0.0 ;
    this->_meanCentred = false ;
  }
  
  
	/*! Assignment operator overload function */
  void operator=(const CMatrixWithMeans &inputMeans )
  { 
		 if (this->pvvd_matrix != NULL) 
			 delete this->pvvd_matrix ;

		 this->pvvd_matrix = (std::vector<std::vector<double> > *) new std::vector<std::vector<double> >() ; 

		 // deep copy of vectors
		 for (size_t t1 = 0 ; t1 < inputMeans.pvvd_matrix->size() ; t1++)
			 this->pvvd_matrix->push_back( ( *(inputMeans.pvvd_matrix))[t1])  ;  // copies the input vectors

     // shallow copy will not assign new vvd structure.
		 // this->pvvd_matrix = inputMeans.pvvd_matrix;
		 //(*(pvvd_matrix))[0].push_back(2.0) ;  // this will add the 2.0 value to the input and the new vectors data


     this->vd_means = inputMeans.vd_means;
     this->_scale  = inputMeans._scale;
     this->_meanCentred  = inputMeans._meanCentred;
  }
  
	std::vector<std::vector<double> > * GetVVDMatrixDataPtr( void ) 
	{
		return this->pvvd_matrix ;
	}

  /*! Calculates the mean of a single vector<double> of the vvd_matrix */
  double CalculateMean(std::vector<double>& input_vect, size_t mean_num )
  {

     //  double s = std::accumulate(input_vect.begin(),input_vect.end(), 0.0);
	  double s = 0.0 ;
	  for (size_t t1 = 0 ; t1 < input_vect.size() ; t1++)
		  s += input_vect[t1] ;

      s = s / (double) input_vect.size() ;
      return s ;
  }
  

  void CalculateMeanAndStore(std::vector<double>& input_vect, size_t mean_num )
  {
      // calculate mean
      double s = CalculateMean(input_vect, mean_num) ;
      
      // Resize storage
      if (this->vd_means.capacity() < this->pvvd_matrix->size() )
      {
        this->vd_means.reserve(this->pvvd_matrix->size());
      }
      // insert mean
      if (mean_num < this->vd_means.capacity() )
      {
          std::vector<double>::iterator d_it;
          d_it = vd_means.begin();
          std::advance(d_it, mean_num) ;
          vd_means.insert ( d_it ,s  );
      }
  }

	/*! 
	 *  Method to copy over CExperimentWithPCAData scores values corresponding to unique batches into the vvd_matrix array and to calculate the mean of the batches.
	 *  This method has an inverse that populates a CExperimentWithPCAData->CExperimentData llpvd array with data from CMatrix vvd array
	 * \param CExperimentWithPCAData The scale factor
	 * \param vd_means_ The vector of means
	 */
	void CreateCMatrixFromCExperimentData(std::vector<std::vector<double> > * pvvd_matrixIN )
  {
		
		if (this->pvvd_matrix != NULL)
			delete this->pvvd_matrix ;

		// assign member structure the input ptr
		this->pvvd_matrix = pvvd_matrixIN ;

    for (size_t batch_int =0 ; batch_int < this->pvvd_matrix->size() ; batch_int++)
    {
      CalculateMeanAndStore( (*(pvvd_matrix))[batch_int], batch_int  ) ; 
    }
  }

  /*! 
	 *  Method to copy over CExperimentWithPCAData scores values corresponding to unique batches into the vvd_matrix array and to calculate the mean of the batches
	 * \param CExperimentWithPCAData The scale factor
	 * \param vd_means_ The vector of means
	 
	void CreateCExperimentDataFromCMatrix(CExperimentWithPCAData * ptr_ExptIN )
  {
		ptr_ExptIN->SetExptDataObjectData(&this->vvd_matrix) ;
  }
	*/

	/*! 
	 * \brief Method to add back a scaled amount of the input vector values to the matrix 'columns'.
	 * \param scale_ The scale factor
	 * \param vd_means_ The input vector of means
	 */
  void ScaleVectors(double scale_, std::vector<double>& vd_means_ )
  {
    this->_scale = scale_ ; 
    // get pointer to first element in input vector
    double * ptr_mean_dbl = &vd_means_[0] ;

    // for each batch present, create a vector containing the samples it holds
    for (size_t i =0 ; i < pvvd_matrix->size() ; i++)  
    {      
      std::vector<double>& v_temp = (*(pvvd_matrix))[i] ; // get the 'i'th vector<double> 
      double * ptr_temp_dbl = &v_temp[0];                 // get a pointer to the first value in the vector<double>
      
      for (size_t j =0 ; j < v_temp.size() ; j++)         // scale each value in the vector<double>
      {
        ptr_temp_dbl[j] += scale_ * ptr_mean_dbl[i] ;
      }
    }
  }
  
  /*! 
   * \brief Counts the sum of the lengths of the component vectors.
   * \return The sum of the sizes of the vectors in pvvd_matrix
   */
  size_t CountValsInVectors( void )
  {
    size_t retcount ;
    retcount = 0 ;
    for (size_t i =0 ; i < pvvd_matrix->size() ; i++)  
    {      
      //std::vector<double>& v_temp = (*(pvvd_matrix))[i] ;
      //retcount += v_temp.size() ;
      retcount += (*(pvvd_matrix))[i].size() ;
    }
    return (retcount) ;
  }
  
  
  /*! 
   * \brief Copies the data in the pvvd_matrix vectors to a single std::vector<double>
   * \return A pointer to the created vector of copied data
   */
  std::vector<double> * ReturnAsVectorOfDoubles( void )
  {
    
    size_t retcount = CountValsInVectors() ;  // get the combined size of all the vectors
    std::vector<double> * ret_vd = (std::vector<double> *) new std::vector<double>( retcount )  ;  // create the return vector
    double * ret_vd_ptr = ret_vd->data() ;
    
    // for each vector present, copy the data to the return vector 
    for (size_t i =0 ; i < pvvd_matrix->size() ; i++)  
    {      
      std::vector<double>& v_temp = (*(pvvd_matrix))[i] ; // get the 'i'th vector<double> 
      double * ptr_temp_dbl = &v_temp[0];                 // get a pointer to the first value in the vector<double>
      
      for (size_t j =0 ; j < v_temp.size() ; j++)         // copy each value in the ith vector<double> to the return vector<double>
      {
        ret_vd_ptr[j] =  ptr_temp_dbl[j] ;
      }
      
      ret_vd_ptr += v_temp.size() ;  // increment the pointer into the return vector
    }
    return (ret_vd) ;
  }
  
  /*! 
   * \brief Copies the data in the pvvd_matrix vectors to a single std::vector<double>
   * \return A pointer to the created vector of copied data
   */
  std::vector<double> * ReturnAsVectorOfDoubles_TEST( void )
  {
    size_t retcount = CountValsInVectors() ;  // get the combined size of all the vectors
    std::vector<double> * ret_vd = (std::vector<double> *) new std::vector<double>( retcount )  ;  // create the return vector
    double * ret_vd_ptr = ret_vd->data() ;
    
    // for each vector present, copy the data to the return vector 
    for (size_t i =0 ; i < pvvd_matrix->size() ; i++)  
    {      
      std::vector<double>& v_temp = (*(pvvd_matrix))[i] ; // get the 'i'th vector<double> 
     //double * ptr_temp_dbl = &v_temp[0];                 // get a pointer to the first value in the vector<double>
      
      for (size_t j =0 ; j < v_temp.size() ; j++)         // copy each value in the ith vector<double> to the return vector<double>
      {
        ret_vd_ptr[j] =  (double) i ; //ptr_temp_dbl[j] ;
      }
      
      ret_vd_ptr += v_temp.size() ;  // increment the pointer into the return vector
    }
    return (ret_vd) ;
  }
  
  

  void PrintMatrixVectors( )
  {
    for (size_t i =0 ; i < pvvd_matrix->size() ; i++)  // for each batch present, create a vector containing the samples it holds
    {      
      std::vector<double>& temp = (*(pvvd_matrix))[i] ;  // get the vector<double> 
      double * ptr_temp_dbl = &temp[0];      
      for (size_t j =0 ; j < temp.size() ; j++)  // for each double in the batch vector
      {
        _PRINTSTD << ptr_temp_dbl[j] << " "  ;
      }
      _PRINTSTD << std::endl ;
    }
  }
  
  void PrintMatrixMeans( )
  {
    //double * ptr_temp_dbl = &this->vd_means[0];
    for (size_t i =0 ; i < this->vd_means.size() ; i++)  // for each batch present, create a vector containing the samples it holds
    {             
        _PRINTSTD << this->vd_means[i] << " "  ;
    }
    _PRINTSTD << std::endl ;
    
    if (this->_meanCentred == true)  
      _PRINTSTD << "This matrix has been mean centred! so printed means are of the original matrix data "  << std::endl ; 
    
  }
  
  void MeanCentre( std::vector<double>& vd_means_ )
  {
     ScaleVectors(-1.0 , vd_means_ ) ;
     this->_meanCentred = true ;
  }




  std::vector<double> ReturnCurrentMeans()
  {
    this->vd_means.clear() ;
    for (size_t batch_int =0 ; batch_int < this->pvvd_matrix->size() ; batch_int++)
    {
      CalculateMeanAndStore( (*(pvvd_matrix))[batch_int], batch_int  ) ; 
    }
    return this->vd_means ; 
  }
  
   std::vector<double> ReturnStoredMeans()
  {
    return this->vd_means ; 
  }
  
} ;  // end class CMatrixWithMeans


/*! Batch probability calculation. Integrates the normal distribution curve to return the confidence (probability) that a batch effect is present.
\brief Returns the confidence (from 0.0 to 1.0) that a batch effect is present.
\note Currently uses an approximation to the cumlative normal proabability function in method norm_cdf_approx(), which has errors of max ~ +-0.0032
*/
class  CCalaculateBatchConfidence
{
	/*! _k_scale_fact is a probability scaling factor which is included to make the 
     P(no batch effect|no batch effect)=1. In the case of no batch effect 
     all the batch means should be equal and their sum is zero as well(mean 
     of each principal component is zero).Then all batch means are zeros and
     the cdf of a normal distribution at 0 is equal to 0.5.Then each of the 
     z products are(0.5)^(b-1) resulting in a sum of products of b*(0.5)^(b-1)
    */
	double _k_scale_fact ; 
	CMapSelectKFromN * _newNKMap ;
	size_t _num_batches ;
	size_t _total_combs ;
	std::vector<double> _z ;
	double _zsop ;  // sum of products
	double _confidence ;

public :

	/*! Constructor */
	CCalaculateBatchConfidence( size_t num_batchesIN, size_t arrayAlignIN ) 
	{
		this->_num_batches = num_batchesIN ;
		this->_k_scale_fact =  pow(2.0, (double)(_num_batches-1)) / (double) (_num_batches) ;

		// The CMapSelectKFromN object lets us access array indecies of a binomial sequence
    // We may want to only create this once and add N,K pairs to the map as different ones are needed ( N,K == 300,7 took about 0.3 seconds to create)
    this->_newNKMap = (CMapSelectKFromN *) new CMapSelectKFromN(std::make_pair(_num_batches,_num_batches-1), arrayAlignIN) ;  // this object lets us access array indecies of a binomial sequence

		this->_total_combs =  this->_newNKMap->nchoosek(_num_batches,_num_batches-1) ;  // this should always be the same as num_batches
		this->_zsop = 0.0 ;
		this->_confidence = -1.0 ;

	}

	~CCalaculateBatchConfidence( void ) 
	{
		delete this->_newNKMap ;
	}


	/*! The equivalent of batchprobability_unbalanced.m from original Matlab code
	*/
	double compute_confidence( std::vector<double> * batch_means, std::vector<double> * sim_means , std::vector<double> * sim_stdevs ) 
	{
		this->_z.reserve(batch_means->size()) ;
		this->_z.assign(batch_means->size(),0.0) ;
		// for (size_t t1 = 0 ; t1 < batch_means->size() ; t1++) this->_z.push_back(0.0) ;
		double * ptr_z = this->_z.data()			;
		double * t_bm  = batch_means->data(); // get pointer to data array of current batch means
		double * t_sm  = sim_means->data()  ; // get pointer to data array of simulated means
		double * t_ss  = sim_stdevs->data()	; // get pointer to data array of simulate standard deviations
		for (size_t t1 = 0 ; t1 < batch_means->size() ; t1++) // for each batch
		{
			// create the NCDF values
			double t_norm_and_scaled  ;
			t_norm_and_scaled = -1.0 * std::abs(t_bm[t1]) ;
			t_norm_and_scaled -=	t_sm[t1] ;
			t_norm_and_scaled /=	t_ss[t1]  ;
			double t_ncdf2  ;
			t_ncdf2 = norm_cdf_approx_closer( t_norm_and_scaled ) ; // return the probability integral
			ptr_z[t1] = t_ncdf2 ;
		}

		double current_prod ;
		this->_zsop = 0.0 ; // set "sum of products" to zero
		std::vector< size_t > v_index(this->_num_batches-1) ; // one for each treatment to be choosen (size() == k)
		// create a set of all combinations of the above _z values (equivalent of Matlab combnk() function)
		for (size_t seq_val = 1 ; seq_val <= this->_total_combs ; seq_val++)
		{            
			// get an array of indexes
			v_index = _newNKMap->GetScan_shortcut_debug(seq_val, (std::string)" compute_confidence() ") ; // requires a sequence that starts at 1 
			current_prod = 1.0 ;
			// select the data from _z via the returned indexes (GetScan returns a binomial distribution sequence)
			for (size_t sum = 0 ; sum < v_index.size() ; sum++)
			{
				current_prod *= _z.at(v_index.at(sum)) ;  // create the product of the values directly
			}
			this->_zsop += current_prod ; // add to sum of products result
		}

		this->_confidence = 1.0 - (this->_k_scale_fact * this->_zsop) ;

		return ( this->_confidence ) ; // return the confidence value
	}

	double sign_ ( double xIN ) 
	{
		return round((xIN / abs(xIN))) ;
	}

	/*! There is a Harman R function that shows the errors between this fast method (from wikipeadea)
	R code: 
		err1 <- norm_cdf_approx_err(seq(-5,5,0.002),0,1)
		plot( seq(-5,5,0.002), err1) 
	This above code shows that the magnitude maximum error is about 0.00031 between the approximate function and the R version of the same function (16 digit accuracy) 
	from : https://en.wikipedia.org/wiki/Normal_distribution 
	*/
	double norm_cdf_approx( double x ) 
	{
			return (0.5 * ( 1 + sign_(x) * sqrt( (1 - exp( (-2 * x * x) / M_PI )) )  ) ) ;
  }



	/*! 
	Ref: Zelen & Severo (1964)  magnitude maximum error is about  |Err(x)| < 7.5x10-8   
	from : https://en.wikipedia.org/wiki/Normal_distribution 
	It seems like the input probability must be +ve to get the accuratcy stated, so we test and correct using the complement of the negated input.
	*/
	double norm_cdf_approx_closer( double x ) 
	{
			
		double phi ;  // standard normal distribution
		double Phi ;  // culmulative standard normal distribution
		double t ;
		const double  b0 = 0.2316419 ;
		const double  b1 = 0.319381530 ;
		const double  b2 = -0.356563782 ;
		const double  b3 = 1.781477937 ;
		const double  b4 = -1.821255978 ;
		const double  b5 = 1.330274429 ;	

		bool x_is_negative = false ;

		x_is_negative = (x < 0) ;

		if (x_is_negative) x = -x ;  // make positive

		t = 1.0 / (1.0 + (b0 * x)) ;
		phi = ( 1.0 / sqrt(2.0 * M_PI)) * exp(-0.5*x*x) ;

    Phi = 1.0 - ((phi * (b1*t))  + (phi * (b2 * t * t))  + (phi * (b3* t * t * t)) + (phi * (b4 * t * t * t * t)) + (phi * (b5* t * t * t * t * t)))  ;
	/* 
		double barray[] = { 0.2316419, 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429}  ;
		double currsum, curr_t ;
		currsum	= 0.0 ;
		curr_t =  t ;
		for (size_t t1 = 1 ; t1 < 6 ; t1++ )
		{			
			currsum += phi * barray[t1]  *curr_t ;
			curr_t *= t ;
		}
		Phi = 1.0 - currsum ;
	*/

		if (x_is_negative) Phi = 1.0 - Phi ;

		return (Phi) ;
		
	}
} ;


/*! \brief Implementation of algorithm to modify batch score values of a single PC/variable to be within defined probablity of having batch effects.
 *  Implements what was present in MATLAB file batchcleaner_unbalancedYO.m
 * 
 */
class  CProcessScoreData 
{
private:
	std::vector<double> _scores_of_PC ;  // this is the input, from the constructor
	std::vector<double> _scores_of_PC_corrected ;  /// this is the output
	std::vector<double> _columnmeans ;  // the original means of the of the batches
	std::vector<double> _batchnmeans_0 ;
	std::vector<double> _batchnmeans_1 ;
	std::vector<double> _batchmeans_m ;
	CMatrixWithMeans * _scores_in_batches ; /*! this is _scores_of_PC data cut into vectors of batches and has the means of those batches in a seperate vector<> */
	CMatrixWithMeans * _zmatrix ; // the mean centered version of _scores_in_batches
	CMatrixWithMeans * _mat_0 ;
	CMatrixWithMeans * _mat_1 ;
	CMatrixWithMeans * _mat_m ;
	double _s0, _s1, _sm ;  // these are the scaling factors
	double _confidence_0, _confidence_1, _confidence_m ;
	double _bestcorrection ;
	double _bestconfidence ;
	
	CCalaculateBatchConfidence * ptr_CalcConfidence ;  
	CSimulateBatchDistribution * ptr_SimBatchDist ;  // this requires the above created CExperimentWithPCAData object in its creator

	/* make default creator private */
	CProcessScoreData()  ;

public: 

	void PrintVectorWithLabel(std::vector<double> * v_IN, std::string LabelIN)
	{
		_PRINTSTD << LabelIN << "\t\t" ;
		for (size_t t1 = 0 ; t1 < v_IN->size() ; t1++)
			_PRINTSTD << v_IN->at(t1) << ",\t"  ;
		_PRINTSTD << std::endl ;
	}
  /*!  
   * \param max_combinationsIN Is the upper limit on the total possible combinations that will be explicitly (comprehensively) computed in 
	 CSimulateBatchDistribution::CalculateXMatrix() and CSimulateBatchDistribution::CalculateMMatrix() methods.
   * \param arrayAlignIN is the alignment of the vector data used by CMapSelectKFromN::GetScan_Vectorised(). Possibly set via -D_ARRAYALIGNEMT= definition in Markevars file
	 * \param confidence_limitIN is user applied limit to restrict the difference between our perturbed batch scores data and the original batch data
   * \param CExperimentWithPCAData Contains the experiment struture and a single PC of data from which to simulate the data distribution function with.
   * \reference Matlab code batchcleaner_unbalancedYO.m - this class represents computations involved a single loop (for a single PC of interest)
   */
  CProcessScoreData(  double confidence_limitIN, CExperimentWithPCAData *  ptr_ExptPCADataIN, size_t max_combinationsIN,  size_t arrayAlignIN, size_t seedIN, bool forceRandomIN) 
  {
		_s0 = 0.0 ;
		_s1 = 1.0 ;
		
	//	_PRINTSTD << "CProcessScoreData(): Entry in to function: " << std::endl ;
		
		ptr_SimBatchDist =  (CSimulateBatchDistribution*) new CSimulateBatchDistribution(arrayAlignIN, ptr_ExptPCADataIN, seedIN )  ;
		ptr_SimBatchDist->FillCombinationsMatrix( max_combinationsIN ) ; // totalnocombs_ptr only needs to be called once
		
		this->ptr_CalcConfidence = (CCalaculateBatchConfidence *) new CCalaculateBatchConfidence(ptr_ExptPCADataIN->GetNumBatches(), arrayAlignIN) ;


		_scores_in_batches = new CMatrixWithMeans() ;
		_scores_in_batches->CreateCMatrixFromCExperimentData( ptr_ExptPCADataIN->CreatePVVD_FromCExperimentData() ) ;
		_columnmeans = _scores_in_batches->ReturnStoredMeans() ; // these are used repeatedly
   
		_zmatrix = new CMatrixWithMeans() ;
		_zmatrix->CreateCMatrixFromCExperimentData( ptr_ExptPCADataIN->CreatePVVD_FromCExperimentData()  ) ;
		_zmatrix->MeanCentre(_zmatrix->vd_means) ;
		_mat_0 = new CMatrixWithMeans() ;
		_mat_1 = new CMatrixWithMeans() ;
		_mat_m = new CMatrixWithMeans() ;

		*_mat_0 = *_zmatrix ;  //  copy mean centred matrix of batches
		*_mat_1 = *_zmatrix ;  //  copy mean centred matrix of batches
    
		// add back a scaled version of the original column (batch) means
		_mat_0->ScaleVectors(_s0, _columnmeans ) ; // multiplies each vector element by the input scale_ value (_s0)
		_mat_1->ScaleVectors(_s1, _columnmeans ) ; // multiplies each vector element by the input scale_ value (_s1)
    
		_batchnmeans_0 = _mat_0->ReturnCurrentMeans() ; // these will be zero (or very close to due to round-off error)
		_batchnmeans_1 = _mat_1->ReturnCurrentMeans() ; // these will be the same as the original scores batch means
   
		ptr_ExptPCADataIN->SetExptDataFromBatchVectors( _mat_0->GetVVDMatrixDataPtr() ) ; // set the data for the simulation
		
		
		ptr_SimBatchDist->SimulateBatches(max_combinationsIN, forceRandomIN) ;
		this->_confidence_0 = this->ptr_CalcConfidence->compute_confidence(&_batchnmeans_0, ptr_SimBatchDist->GetPtrMeans(), ptr_SimBatchDist->GetPtrStdDevs() ) ;

#if (_PRINT==1)
		PrintVectorWithLabel(&_batchnmeans_0, "batch means 0:") ;
		PrintVectorWithLabel(ptr_SimBatchDist->GetPtrMeans(), "sim means   :") ;
		PrintVectorWithLabel(ptr_SimBatchDist->GetPtrStdDevs(), "sim stdevs :") ;
#endif

		ptr_ExptPCADataIN->SetExptDataFromBatchVectors( _mat_1->GetVVDMatrixDataPtr() ) ; // set the data for the simulation
		ptr_SimBatchDist->SimulateBatches(max_combinationsIN, forceRandomIN) ;
		this->_confidence_1 = this->ptr_CalcConfidence->compute_confidence(&_batchnmeans_1, ptr_SimBatchDist->GetPtrMeans(), ptr_SimBatchDist->GetPtrStdDevs() ) ;
#if (_PRINT==1)
		PrintVectorWithLabel(&_batchnmeans_1, "batch means 1:") ;
		PrintVectorWithLabel(ptr_SimBatchDist->GetPtrMeans(), "sim means   :") ;
		PrintVectorWithLabel(ptr_SimBatchDist->GetPtrStdDevs(), "sim stdevs :") ;
		//_PRINTSTD << "confidence 0 = " << this->_confidence_0 << std::endl; _PRINTSTD  << "confidence 1 = " << this->_confidence_1 << std::endl ;
#endif
		TOTAL_ITERATIONS +=2;
#if (_PRINT==1)
		_PRINTSTD << "ITERATION:\t" << TOTAL_ITERATIONS << ",\t confidence 0 = \t" << this->_confidence_0   << " \tconfidence 1 = \t" << this->_confidence_1 << std::endl ;
#endif
    size_t count ;
		count = 0 ;
		if ((_confidence_1-confidence_limitIN) > 0 )   // To check if shrinking is needed or not
		{
				while ((_s1-_s0) > 0.01)                 // Go for shrinking till confidence gets to its limit 
				{
						count++ ;
						_sm= (_s1+_s0)/2.0 ;

					*_mat_m = *_zmatrix ;
					 _mat_m->ScaleVectors(_sm, _columnmeans ) ;  // Shrink mean pca score matrix
					 _batchmeans_m = _mat_m->ReturnCurrentMeans() ;
					 ptr_ExptPCADataIN->SetExptDataFromBatchVectors( _mat_m->GetVVDMatrixDataPtr() ) ; // set the data for the simulation
					 ptr_SimBatchDist->SimulateBatches(max_combinationsIN, forceRandomIN) ; // simulate the new data

					_confidence_m = this->ptr_CalcConfidence->compute_confidence(&_batchmeans_m, ptr_SimBatchDist->GetPtrMeans(), ptr_SimBatchDist->GetPtrStdDevs() ) ; 
#if (_PRINT==1)
					PrintVectorWithLabel(&_batchmeans_m, "batch means m:") ;
					PrintVectorWithLabel(ptr_SimBatchDist->GetPtrMeans(), "sim means   :") ;
					PrintVectorWithLabel(ptr_SimBatchDist->GetPtrStdDevs(), "sim stdevs :") ;

#endif
					if (_confidence_m < confidence_limitIN)
					{
						_s0=_sm;
						//*_mat_0 =*_mat_m ;
						_confidence_0=_confidence_m;
					}
					else
					{
						_s1=_sm;
						//*_mat_1 =*_mat_m ;
						_confidence_1=_confidence_m;
					}

					TOTAL_ITERATIONS++ ;
#if (_PRINT==1)
					_PRINTSTD << "ITERATION:\t" << TOTAL_ITERATIONS << ",\t confidence 0 = \t" << this->_confidence_0  << "\tconfidence 1 = \t" << this->_confidence_1 << std::endl ;
#endif
				} // end of while loop

				if ((_confidence_1-confidence_limitIN) < (confidence_limitIN-_confidence_0))
				{
						_bestcorrection=floor(100.0*_s1)/100.0;
						_bestconfidence=_confidence_1;
				}
				else
				{
						_bestcorrection=ceil(100.0*_s0)/100.0;
						_bestconfidence=_confidence_0;
				}

		}
		else
		{
			_bestcorrection=1.0;
			_bestconfidence=_confidence_1;
		}

		
  }

	double GetBestCorrection( void ) 
	{
		return this->_bestcorrection ;
	}

  double GetBestCConfidence( void ) 
	{
		return this->_bestconfidence ;
	}


  ~CProcessScoreData()
  {
      delete ptr_SimBatchDist ;
			delete ptr_CalcConfidence ;

			delete _scores_in_batches ; // this is _scores_of_PC data cut into vectors of batches and has the means of those batches in a seperate vector<>
      delete _zmatrix ; // the mean centered version of _scores_in_batches
      delete _mat_0 ;
      delete _mat_1 ;
  }
  
  
} ;


#if defined (_MSC_VER)

/***********************************************************************************
 * The following code is used to compile the project source as a stand alone executable
 * There is a Makefile in the Harman/src directory that can be used for compilation, 
 * however it must be renamed for when "R CMD Install" is run as it overrides the Makevars 
 * and Makevars.win files.
 **************************************************************************************/
 
void PrintUsage(std::string progName) 
{
	  _PRINTSTD << "Usage: " << progName << " -d [datapath] -t [scores filename] -b [batch info filename] -f [0|1]] -v [0|1]] -s [samples (== rows)] -p [PCs == columns]" << std::endl ;
		_PRINTSTD <<  " -d [datapath]   Path to the data input files  " << std::endl ;
		_PRINTSTD <<  " -t [scores filename] 	 " << std::endl ;
		_PRINTSTD <<  " -b [batch info filename] 	 " << std::endl ;
		_PRINTSTD <<  " -s [samples]   	Number of samples in the data  " << std::endl ;
		_PRINTSTD <<  " -p [PCs] 				Number of PCs in the dataset " << std::endl ;
		_PRINTSTD <<  " -f [0 | 1]] 				Force random selection code to compute mean " << std::endl ;
		_PRINTSTD <<  " -n select k from n	 " << std::endl ;
		_PRINTSTD <<  " -k select k from n	 " << std::endl ;
		//_PRINTSTD <<  " -o [0 | 1] if 1 then use OpenCL select K form N code " << std::endl ;
		_PRINTSTD << std::endl ;
}


/*! \brief Struct used for to store inputs from command line arguments using the getopt.h library
 Struct used for to store inputs from command line arguments using the getopt.h library*/
struct arg_list {
	double _limit ;										/*! [-l] The maximum error limit for the approximation */
	std::string _datapath  ;					/*! [-d] Directory where the input files are stored */
	std::string _scores_filename ;		/*! [-t] Input data file containing the scores. Dimension = [samples x PCs] == [rows x columns] */
	std::string _batchinfo_filename ; /*! [-b] Description of what batch has what treatment. Comprises of 2 rows: 1 for batch designation and the other for treatment designation */
	size_t _samples ;									/*! [-s] The number of samples in the experiment == scores rows (default = 24)  */
  size_t _PCs ;											/*! [-p] The number of principal compoenents (PC's) in the input scores file == scores columns (default = 23) */
	bool _forceRandom  ;							/*! [-f] Force the random selection of data to get the mean */
	bool _printbool  ;								/*! [-v] verbose output to screen */
	bool USEOPENCL_SELECT_ARG ;  /*! Used in GPU version */
	size_t N ; /*! Used to test select K from N codes */
	size_t K ; /*! Used to test select K from N codes */
	size_t numrepeatsKfromN ; /*! Used to test selct K from N codes */
};


//struct arg_list get_opencl_args(int argc, wchar_t * const* argv)
struct arg_list get_opencl_args(int argc, char* argv[])
{
	int opt;
	struct arg_list ret_list ;

	// set the defaults
	ret_list._limit = 0.95 ;
	ret_list._printbool = false ;
	ret_list._datapath.assign(_DATAPATH) ;
	ret_list._batchinfo_filename.assign("bladder_expt_batch.bin") ;  // batch_balanced_npm
  ret_list._scores_filename.assign("bladderdata_scores.bin") ;
	ret_list._samples =  57 ;
	ret_list._PCs = 56 ;
	ret_list._forceRandom = false ;
	ret_list.USEOPENCL_SELECT_ARG = false ;
	ret_list.N = 20 ;
	ret_list.K = 4 ;
	ret_list._printbool = true ;
	ret_list.numrepeatsKfromN = 100000 ;
	// "Usage: a.out -d [datapath] -t [scores input data file ] -b [batch and treatment info file] -s [number of samples] -p [number of PCs == columns] -v [verbose]"
	// Parse command-line arguments
	 const char* optString = "d:l:s:p:t:b:n:k:o:f:r:vh";
	//const char* optString = "d:s:p:h";
	opt = getopt(argc, argv, optString);
	/*if (opt == -1) {
		PrintUsage(argv[0]) ;
		_PRINTERROR <<  "User supplied arguments: " << std::endl ;
		for (int t1 = 1; t1 < argc; t1++)
			_PRINTERROR <<  argv[t1] << ", " ;
		exit(EXIT_FAILURE) ;
	}*/
	while (opt != -1)
	{
		switch (opt) {
		case ('?'):  // print usage message
			PrintUsage(argv[0]) ;
			exit(EXIT_FAILURE) ;
			break ;
		case ('h'):  // print usage message
			PrintUsage(argv[0]) ;
			exit(0) ;
			break ;
		case ('d'):  // get filesystem path to data string
			//ret_list.plat_str = optarg ;
			ret_list._datapath.assign(optarg);
			break ;
		case ('l'):  // get filesystem path to data string
			//ret_list.plat_str = optarg ;
			ret_list._limit = 0.95 ;
			break ;
		case ('t'):  // get scores file data string
			//ret_list.plat_str = optarg ;
			ret_list._scores_filename.assign(optarg);
			break ;
		case ('b'):  // get batch information file string
			//ret_list.plat_str = optarg ;
			ret_list._batchinfo_filename.assign(optarg);
			break ;
		case ('s'):  // get number of samples
			ret_list._samples = strtoll(optarg, NULL, 10) ;
			break ;
		case ('v'):  // verbose output where available
			if (strtoll(optarg, NULL, 10) == 1)
				ret_list._printbool = true ;
			else 
				ret_list._printbool = false ;
			break ;
		case ('p'):// get number of PCs
			ret_list._PCs = strtoll(optarg, NULL, 10) ;
			break ;
		case ('f'):// force random
			if (strtoll(optarg, NULL, 10) == 1)
				ret_list._forceRandom = true ;
			else 
				ret_list._forceRandom = false ;
			break ;
		case ('n'):// select K from 'N' 
			ret_list.N = strtoll(optarg, NULL, 10) ;
			break ;
		case ('k'):// select 'K' from N
			ret_list.K = strtoll(optarg, NULL, 10) ;
			break ;
		case ('r'):// the number of repeates for testing selectKfromN code
			ret_list.numrepeatsKfromN = strtoll(optarg, NULL, 10) ;
			break ;
		case('o') :
			if (strtoll(optarg, NULL, 10) == 1)
				ret_list.USEOPENCL_SELECT_ARG = true ;
		
		}
		opt = getopt(argc, argv, optString) ;
	}
	return (ret_list) ;
}


// used for debugging via C++ debugger
/*! main() entry point to C/C++ version of the code.
 * 
 *  This entry point is used if _MSC_VER is defined (which we can use as a hack on Linux also via -D_MSC_VER=1 on the compilation command line)
 */
// g++ -m64 -I"F:/R/R-3.0.2/include" -DNDEBUG -IF:/R/library/Rcpp/include -I"F:/R/library/Rcpp/include" -I"d:/RCompile/CRANpkg/extralibs64/local/include"  -O0 -g -D_MSC_VER=1 -Wall  -mtune=core2 -c Harman.cpp -o Harman.exe -Ld:/RCompile/CRANpkg/extralibs64/local/lib/x64 -Ld:/RCompile/CRANpkg/extralibs64/local/lib -LF:/R/R-3.0.2/bin/x64 -lR

int main(int argc, char* argv[])
{
	struct arg_list main_args ;
	main_args = get_opencl_args(argc, argv) ;
	//main_args._batchinfo_filename.assign("batch_balanced_npm.bin") ;  // batch_unbalanced_npm

  std::vector<double>  scores_std ; //!< points to the same data as in R (use clone<>() to make a copy)  
  std::vector<size_t>  tr_std ;
  std::vector<size_t>  ba_std ;
	tr_std.reserve(main_args._samples) ;
	ba_std.reserve(main_args._samples) ;
	scores_std.reserve(main_args._samples * main_args._PCs) ;

	// select file from disk
	std::string istream_filename ;
	istream_filename = main_args._datapath + main_args._batchinfo_filename  ;
	_PRINTSTD << "opening file: " << istream_filename << std::endl ;
	std::ifstream filein(istream_filename.c_str(), std::ios::binary) ;
	int tint ;
	if (filein.is_open()) {
		while (filein) {
			// R saves matrix data in column major order so the batch and treatment info will be 'interlaced'
			// N.B. the scores data was transposed before saving
			for (int t1 = 0 ; (t1 < (main_args._samples)) && (filein) ; t1++)
			{
				filein.read((char*)&tint, sizeof(int)) ;
   // Uncomment if we want a 'balanced' version of the NPT data set
	//			if (t1 == 6) tint += 1 ;
	//			if (t1 == 7) tint += 1 ; 
				tr_std.push_back(tint);
	//			filein.read((char*)&tint, sizeof(int)) ;
	//			ba_std.push_back(tint);
			}
			for (int t1 = 0 ; (t1 < (main_args._samples)) && (filein) ; t1++)
			{
			//	filein.read((char*)&tint, sizeof(int)) ;
			//		tr_std.push_back(tint);
				filein.read((char*)&tint, sizeof(int)) ;
				ba_std.push_back(tint);
			}
			break ;
		}
	} 
	else {
		_PRINTERROR << "Error main(): filename could not be opened " << istream_filename << std::endl ;
	}
	filein.close() ;


  // Test some of the random selection code:
	_PRINTSTD << "C++ version: " << __cplusplus << std::endl ;
	_PRINTSTD << "RAND_MAX : " << RAND_MAX << std::endl ;


	_PRINTSTD << "BATCHES : \t" ;
	for (int t1 = 0 ; t1 < main_args._samples ; t1++)
	{
		_PRINTSTD  << ba_std[t1]  << ",\t";
	}
	_PRINTSTD << std::endl ;
	_PRINTSTD << "TREATMENTS :\t" ;
	for (int t1 = 0 ; t1 < main_args._samples ; t1++)
	{
		_PRINTSTD  << tr_std[t1]  << ",\t";
	}
	_PRINTSTD << std::endl ;

	// read the scores data:
	istream_filename = main_args._datapath + main_args._scores_filename  ;
	_PRINTSTD << "opening file: " << istream_filename << std::endl ;
	std::ifstream filein2(istream_filename.c_str(), std::ios::binary) ;
	double t_dbl ;
	if (filein2.is_open()) {
		while (filein2) {
			for (int t1 = 0 ; t1 < main_args._samples * main_args._PCs ; t1++)
			{
				filein2.read((char*)&t_dbl, sizeof(double)) ;
				scores_std.push_back(t_dbl);
			}
		}
	} 
	else {
		_PRINTERROR << "Error main(): filename could not be opened " << istream_filename << std::endl ;
	}
	filein2.close() ;
	

	TOTAL_ITERATIONS = 0 ;
	// create the data structure for the PCA data
	// pcascores pointer to class that holds a copy of the full Scores input data and  its dimensions and structure (row are samples, columns are variables or vice-versa) 
	CPCAScoresArray * pcascores								= (CPCAScoresArray*) new CPCAScoresArray(scores_std, main_args._samples,  main_args._PCs, SamplesAsCols) ;  // creates the pcascores object which stores all the PCA scores for an experiment      
	// Create the Experiment structure object
	CExperimentStructure *  ptr_Expt					= (CExperimentStructure *) new CExperimentStructure( tr_std, ba_std, main_args._samples) ;
	

	size_t maxpc ;
	maxpc = main_args._PCs ;

 ////////////////////// // 
 // This is the new version of the loop:
  size_t calc_seed ;
  // Progress p(nb, true);
  std::vector<double> confidence(main_args._PCs) ;
  std::vector<double> correction(main_args._PCs) ;
  double * confidence_ptr = confidence.data() ;
  double * correction_ptr = correction.data() ;
  
  // Available through library(RcppProgress):
  //  if ( progmon == NULL )
  //     progmon = (Progress * ) new Progress(PCs, false);  // set to true if we want to display a progress bar
  
  bool abort_bool  ;
  abort_bool  = false ;
  
//  #pragma omp parallel for  shared(ptr_Expt, pcascores, abort_bool) ordered schedule(static,1)
  for (size_t numPCs = 0 ; numPCs < main_args._PCs ; numPCs++ )
  {
    bool resbool ;
    resbool =  abort_bool == false ;
    if (resbool )
    {
      
      //int tid = omp_get_thread_num() ; // requires #include <omp.h>
      std::vector<double> * t_PCAScores = pcascores->GetPCData(numPCs) ;
      
      CMatrixWithMeans *    corrected_scores_as_vvd ;
      std::vector<double>   corrected_scores_columnmeans ;
      std::vector<double> * corrected_scores_as_vd ;
      
      CExperimentWithPCAData *  ptr_ExptPCAData = (CExperimentWithPCAData *) new CExperimentWithPCAData(ptr_Expt, t_PCAScores , numPCs) ;
      
      // Get a copy here before things are modified
      corrected_scores_as_vvd = new CMatrixWithMeans() ;
      corrected_scores_as_vvd->CreateCMatrixFromCExperimentData( ptr_ExptPCAData->CreatePVVD_FromCExperimentData() ) ;
      corrected_scores_columnmeans = corrected_scores_as_vvd->ReturnStoredMeans() ; 
      // Mean centre the matrix, as we will add back the corrected mean value
      corrected_scores_as_vvd->MeanCentre(corrected_scores_as_vvd->vd_means) ;
      
      // Use calc_seed to make sure the random seed is different for each iteration (or if randseed is 0, then the system clock is used - which may also casue issues due to low time resolution)
      calc_seed = (numPCs+1) * 42 ; 
      // determine the correction factor if needed 
      CProcessScoreData * ptr_ProcsessScore = (CProcessScoreData *) new CProcessScoreData(0.95, ptr_ExptPCAData, main_args.numrepeatsKfromN, 4, calc_seed, main_args._forceRandom) ;
      
      // Add back the corrected mean
      corrected_scores_as_vvd->ScaleVectors(ptr_ProcsessScore->GetBestCorrection(), corrected_scores_columnmeans ) ;
      //corrected_scores_as_vd = corrected_scores_as_vvd->ReturnAsVectorOfDoubles() ; 
      corrected_scores_as_vd = corrected_scores_as_vvd->ReturnAsVectorOfDoubles_TEST() ; 
      
      // place corrected scores data back into the scores matrix
      //pcascores->SetPCData_TEST(numPCs, corrected_scores_as_vd) ;
      pcascores->SetPCData(numPCs, corrected_scores_as_vd) ;
      
      confidence_ptr[numPCs] = ptr_ProcsessScore->GetBestCConfidence() ;
      correction_ptr[numPCs] = ptr_ProcsessScore->GetBestCorrection() ; 
      
      if (main_args._printbool) {
        _PRINTSTD << "PC :\t" << numPCs+1 << "/" << main_args._PCs << "\tconfidence: \t" << ptr_ProcsessScore->GetBestCConfidence() << " \tcorrection: \t" << ptr_ProcsessScore->GetBestCorrection() << std::endl ;
      }
      
      delete t_PCAScores ;
      delete ptr_ExptPCAData ;
      delete ptr_ProcsessScore ;
      delete corrected_scores_as_vvd ;
      delete corrected_scores_as_vd ;
    }  
    
    
  } // end of openmp for loop

  ////////////////////// // 
  // end of new version of the loop
  ////////////////////// // 

  
	_PRINTSTD << "TOTAL_ITERATIONS :" << TOTAL_ITERATIONS << std::endl ;


  return 0 ;
  //return Rcpp::wrap( 0 );
  
} // main()
 

#endif 



#if !defined (_MSC_VER)



/***********************************************************************************
 * The following code is used as the interface to the R functions that are stored in the 
 * Harman/R directory in file TestHarmanFunctions.R .
 * test_randomSelctFromIntegerVect <- function( inputDataIN,  numselectIN,  numrepsIN,  seedIN)
 ***********************************************************************************/  




// Current R interface:


/**
 * \brief The main entry point for Harman functionality. Inputs the PCA scores and the experiment structure (treatments and batches) and returns a batch corrected version of the PCA scores matrix
 *   
 * \details Input of PCA scores and the experiment structure (treatments and batches) and returns a batch corrected version of the PCA scores matrix
 * 
 * \param [in] SEXP pc_scores_rIN    2D NumericMatrix of scores data. rows = samples, cols = PC scores
 * \param [in] SEXP group_rIN        The structure of the experiment, consisting of batch numbers and treatment numbers
 *                                forming 2 rows or columns (HarmanMain works out which). Each entry for a sample
 *                                describes what batch it came from and what treatment it was given.
 *                                Has to be integer formated data.
 * \param [in] SEXP limit_rIN       A double precsion value indicating the limit of confidence in which to stop removing a batch effect
 * \param [in] SEXP numrepeats_rIN  The number of repeats in which to run the simulated batch mean distribution estimator. Probably should be greater than 100,000.
 * \param [in] SEXP randseed_rIN    Random seed to pass to the random number generator (0 for use default from system time)
 * \param [in] SEXP forceRandomIN   Force the simulated mean code to use random selection of scores to create the simulated batch mean (rather than full explicit calculation from all permutations).
 * \return    SEXP  R list:       scores.corrected  = harman_res_list["corrected_scores"]
 *                                    correction    = harman_res_list["correction"]
 *                                    confidence    = harman_res_list["confidence"]
 *                                    
 * test with:
Windows:
  library(Harman)
  batches <- read.table("F:\\projects\\harman_git\\data\\harman\\Rob_Dunne_info.csv", header=TRUE, sep=",", row.names="Sample")
  datamatrix <- read.table("F:\\projects\\harman_git\\data\\harman\\Rob_Dunne_data_original.csv", header=TRUE, sep=",", row.names="probeID")
 
 batches <- read.table("F:\\projects\\harman_git\\data\\harman\\NPM_info.csv", header=TRUE, sep=",", row.names="Sample")
 datamatrix <- read.table("F:\\projects\\harman_git\\data\\harman\\NPM_data.csv", header=TRUE, sep=",", row.names="probeID")
 
 expt_<- as.factor(batches$Treatment)
 batch_ <- as.factor(batches$Batch)
 system.time(res  <- harman(as.matrix(datamatrix), expt=expt_, batch=batch_, strict=FALSE, printInfo=T) )
 
#system.time(res <- HarmanRun(as.matrix(datamatrix), batches, 0.95, 100000, 42, F, T, strict=FALSE) )
 
 
Linux:
                      
 batches <- read.table("/home/bow355/ACP_HARMAN2015/harman_git/data/Rob_Dunne_info.csv", header=TRUE, sep=",", row.names="Sample")
 datamatrix <- read.table("/home/bow355/ACP_HARMAN2015/harman_git/data/Rob_Dunne_data_original.csv", header=TRUE, sep=",", row.names="probeID")
 expt_<- batches$Treatment
 batch_ <- batches$Batch
 system.time(res  <- harman(as.matrix(datamatrix), expt=expt_, batch=batch_) )
# res <- HarmanRun(as.matrix(datamatrix), batches, 0.95, 100000, 42, F, T, strict=TRUE) 
 
 
 * \useDynLib Harman
 * \importFrom Rcpp sourceCpp
 * 
 */

SEXP HarmanMain(SEXP pc_scores_rIN, SEXP group_rIN, SEXP limit_rIN, SEXP numrepeats_rIN, SEXP randseed_rIN, SEXP forceRandom_rIN, SEXP printInfoIN)
{
  Rcpp::NumericMatrix scores(pc_scores_rIN) ; //!< points to the same data as in R (use clone<>() to make a copy)
  Rcpp::NumericMatrix group(group_rIN) ;   //!< points to the same data as in R (use clone<>() to make a copy)

  double limit                = Rcpp::as<double>(limit_rIN) ; 
  unsigned int numrepeats     = Rcpp::as<unsigned int>(numrepeats_rIN) ;
  unsigned int randseed       = Rcpp::as<unsigned int>(randseed_rIN) ; 
  bool forceRandom            = Rcpp::as<bool>(forceRandom_rIN) ;
  bool printInfo              = Rcpp::as<bool>(printInfoIN) ;
  
  int samples ;
  int PCs ;
  
  // TODO: N.B. hese swap depending on where data arrives from - .bin date has PCs = ncol  and .csv data has PCs = nrow
  // may need to specify SamplesAsRows/SamplesAsCols to select
        samples               = scores.nrow() ;
        PCs                   = scores.ncol() ;
  
  if (printInfo)
  {
    Rcpp::Rcout << "HarmanMain() called with scores matrix dimension: " <<  std::endl ;
    Rcpp::Rcout << "HarmanMain()  Samples \t "  << samples <<  std::endl ;
    Rcpp::Rcout << "HarmanMain()  PCs     \t "  << PCs << std::endl ;
    Rcpp::Rcout << "HarmanMain() called with limit: \t\t\t" << limit << std::endl ;
    Rcpp::Rcout << "HarmanMain() called with numrepeats:\t\t\t" << numrepeats << std::endl ;
    Rcpp::Rcout << "HarmanMain() called with randseed: \t\t\t" << randseed << std::endl ;
    Rcpp::Rcout << "HarmanMain() called with forceRandom: \t\t\t" << forceRandom << std::endl ;
  }
  
  
  Rcpp::NumericVector tr ;
  Rcpp::NumericVector ba ;
  
  // Determine the layout of the input batch information. 
  // group should have either 2 rows or 2 columns and the other index should be the size of the number of samples present.
  if (group.nrow() == 2)
  {
    tr = Rcpp::clone<Rcpp::NumericVector>(group( 0, Rcpp::_)) ;  // treatment data in a row is assumed
    ba = Rcpp::clone<Rcpp::NumericVector>(group( 1, Rcpp::_)) ;  // batch data in a row is assumed
  }
  else
    if (group.ncol() == 2)
    {
      tr = Rcpp::clone<Rcpp::NumericVector>(group(  Rcpp::_,0)) ;  //  treatment data in a  col is assumed
      ba = Rcpp::clone<Rcpp::NumericVector>(group(  Rcpp::_,1)) ;  //  batch data in a  col is assumed
    }
    else
    {
      Rcpp::Rcout <<  "Error: HarmanMain() Input \'group\' should have either 2 rows or 2 columns of batch information data"  <<  std::endl ;
      Rcpp::Rcout << " group.ncol() = " << group.ncol() << " and group.nrow() " << group.nrow() << std::endl ; 
      Rcpp::Rcout << std::flush ;
    }
    
    
    std::vector<unsigned int>    tr_std_tmp     = Rcpp::as<std::vector<unsigned int> >(tr) ;
    std::vector<unsigned int>    ba_std_tmp     = Rcpp::as<std::vector<unsigned int> >(ba) ; 
    std::vector<size_t>    tr_std(tr_std_tmp.size()) ;
    tr_std.assign(&tr_std_tmp[0], &tr_std_tmp[0]+tr_std_tmp.size());
    std::vector<size_t>    ba_std(ba_std_tmp.size()) ; 
    ba_std.assign(&ba_std_tmp[0], &ba_std_tmp[0]+ba_std_tmp.size());

    // I think this takes a copy of the scores data
    std::vector<double> scores_std = Rcpp::as<std::vector<double> >(scores) ; // this may be transposed wrongly - so set SamplesAsRows == SamplesAsCols ?
    // create the data structure for the PCA data
    // pcascores pointer to class that holds a copy of the full Scores input data and  its dimensions and structure (row are samples, columns SamplesAsCols  are variables or vice-versa)
    CPCAScoresArray * pcascores							= (CPCAScoresArray*) new CPCAScoresArray(scores_std, samples,  PCs, SamplesAsCols) ;  // creates the pcascores object which stores all the PCA scores for an experiment
    // Create the Experiment structure object
    CExperimentStructure *  ptr_Expt					= (CExperimentStructure *) new CExperimentStructure( tr_std, ba_std, samples) ;
  
   
   
    int calc_seed ;
   // Progress p(nb, true);
   std::vector<double> confidence(PCs) ;
   std::vector<double> correction(PCs) ;
   double * confidence_ptr = confidence.data() ;
   double * correction_ptr = correction.data() ;
   

   bool abort_bool  ;
   abort_bool  = false ;
   
#if (_USE_RCPP==1)
   GetRNGstate();
#endif
   
   // critical sections around the random (_RANDFUNCTION_) calls in file CSelectRandom.h 
   // casue inefficiencies in this parallel version using R's random number generator
#pragma omp parallel for  shared(ptr_Expt, pcascores, abort_bool) ordered schedule(static,1)
    for (int numPCs = 0 ; numPCs < PCs ; numPCs++ )
    {

        bool resbool ;
        #pragma omp critical 
        {
            resbool =  abort_bool == false ;
        }
      
        if (resbool )
        {
            std::vector<double> * t_PCAScores = pcascores->GetPCData(numPCs) ;

        /* // For testing
          std::vector<double> * t_PCAScores =  (std::vector<double> *) new std::vector<double>(pcascores->_numSamples) ;
          for (int t1 = 0 ; t1 < pcascores->_numSamples; t1++ ) {
            t_PCAScores->at(t1) = t1 ;
          }*/
          
          
            CMatrixWithMeans *    corrected_scores_as_vvd ;
            std::vector<double>   corrected_scores_columnmeans ;
            std::vector<double> * corrected_scores_as_vd ;
            
            CExperimentWithPCAData *  ptr_ExptPCAData = (CExperimentWithPCAData *) new CExperimentWithPCAData(ptr_Expt, t_PCAScores , numPCs) ;
            
            // Get a copy here before things are modified
            corrected_scores_as_vvd = new CMatrixWithMeans() ;
            corrected_scores_as_vvd->CreateCMatrixFromCExperimentData( ptr_ExptPCAData->CreatePVVD_FromCExperimentData() ) ;  // this has a reverse below using SetExptDataFromBatchVectors()
            corrected_scores_columnmeans = corrected_scores_as_vvd->ReturnStoredMeans() ; 
            // Mean centre the matrix, as we will add back the corrected mean value
            corrected_scores_as_vvd->MeanCentre(corrected_scores_as_vvd->vd_means) ;
            
            // Use calc_seed to make sure the random seed is different for each iteration (or if randseed is 0, then the system clock is used - which may also casue issues due to low time resolution)
            calc_seed = (numPCs+1) * randseed ; 
            // determine the correction factor if needed 
            
            // Original- put back when code tested
            CProcessScoreData * ptr_ProcsessScore = (CProcessScoreData *) new CProcessScoreData(limit, ptr_ExptPCAData, numrepeats, 4, calc_seed, forceRandom) ;
           
            // Add back the corrected mean
            corrected_scores_as_vvd->ScaleVectors(ptr_ProcsessScore->GetBestCorrection(), corrected_scores_columnmeans ) ;
            ptr_ExptPCAData->SetExptDataFromBatchVectors( corrected_scores_as_vvd->GetVVDMatrixDataPtr() ) ;
            corrected_scores_as_vd =  ptr_ExptPCAData->Return_vector_from_T() ;
            // place corrected scores data back into the scores matrix
            pcascores->SetPCData(numPCs, corrected_scores_as_vd) ;
            
           
            ///// Testing code version://///
            /*
            corrected_scores_as_vvd->ScaleVectors(1.0, corrected_scores_columnmeans ) ;
            ptr_ExptPCAData->SetExptDataFromBatchVectors( corrected_scores_as_vvd->GetVVDMatrixDataPtr() ) ;  // this is the inverse of CreateCMatrixFromCExperimentData()
            // corrected_scores_as_vd = corrected_scores_as_vvd->ReturnAsVectorOfDoubles_TEST() ; // TEST version returns the batch allocated variables
            corrected_scores_as_vd =  ptr_ExptPCAData->Return_vector_from_T() ;
          //  corrected_scores_as_vd = corrected_scores_as_vvd->ReturnAsVectorOfDoubles() ;
            pcascores->SetPCData(numPCs, corrected_scores_as_vd) ;
            */
            
            confidence_ptr[numPCs] = ptr_ProcsessScore->GetBestCConfidence() ;
            correction_ptr[numPCs] = ptr_ProcsessScore->GetBestCorrection() ; 
           
           if (printInfo) {
            #pragma omp ordered
             Rcpp::Rcout << "PC :\t" << numPCs+1 << "/" << PCs << "\tconfidence: \t" << ptr_ProcsessScore->GetBestCConfidence() << " \tcorrection: \t" << ptr_ProcsessScore->GetBestCorrection() << std::endl ;
           }
           
            delete t_PCAScores ;
            delete ptr_ExptPCAData ;
            delete ptr_ProcsessScore ;
            delete corrected_scores_as_vvd ;
            delete corrected_scores_as_vd ;
        }  
        
        // Rcpp::checkUserInterrupt() ;
        // Temporarily comment out the next omp master pragma to make the code compile.
        
       // #pragma omp master
#ifdef _OPENMP
       if ( omp_get_thread_num()==0 ) // requires #include <omp.h>
#endif
        {
            #pragma omp critical
            {
                    abort_bool = true ;
                #pragma omp flush (abort_bool)
                    Rcpp::checkUserInterrupt() ;
                    abort_bool = false ; // if interupted then we will not get to here via the master thread
                #pragma omp flush (abort_bool)        
            
            } // end critical
        }
         
 
  } // end of openmp for loop
    
#if (_USE_RCPP==1)
    	PutRNGstate();
#endif
    Rcpp::NumericVector corrected_scores(pcascores->_scores.begin(), pcascores->_scores.end()) ;
  /* 
    // Copy stl results into Rcpp vectors:
    Rcpp::NumericVector confidence_rcpp(confidence.begin(), confidence.end() ) ; // = Rcpp::NumericVector::create( PCs ) ;
    Rcpp::NumericVector correction_rcpp(correction.begin(), correction.end() ) ; // = Rcpp::NumericVector::create( PCs ) ;
    Rcpp::NumericVector corrected_scores(pcascores->_scores.begin(), pcascores->_scores.end()) ;
 

    // N.B. SamplesAsCols is the current (C/C++) data format
    double * scores_ptr = pcascores->ReturnPointerToVectorData() ;
    Rcpp::NumericMatrix corrected_scores(PCs, samples) ;  // As R is back to front (Fortran) major. we will swap the R side matrix to have "samples as rows"
    for( size_t pcnum = 0 ; pcnum < PCs ; pcnum++ ){
      for( size_t sample_num=0; sample_num<samples; sample_num++){
        corrected_scores(pcnum, sample_num ) = scores_ptr[pcnum * PCs +  sample_num ] ;
      }
    }
 */
  
  delete pcascores ;
  delete ptr_Expt ;
  
  if (abort_bool == false) {
    return Rcpp::wrap(Rcpp::List::create(
                                          Rcpp::Named("corrected_scores") = corrected_scores,  // pcascores->_scores,  // corrected_scores,
                                          Rcpp::Named("confidence_vector") = confidence,
                                          Rcpp::Named("correction_vector") = correction 
                                        )
                     );
  }
  else
  {
    return Rcpp::wrap(Rcpp::List::create(
                                          Rcpp::Named("corrected_scores") = NULL,
                                          Rcpp::Named("confidence_vector") = NULL,
                                          Rcpp::Named("correction_vector") = NULL
                                        )
                      );
  }
}

  

#endif // #if !defined (_MSC_VER)
