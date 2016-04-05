/**
 * \file Harman/src/CSelectRandom.h 
 * \date 12/11/2015
 * \author   Josh Bowden, CSIRO
 * \details Harman eResearch Collaboration project ERRFP-263 (https://jira.csiro.au/browse/ERRFP-263)
 * Stash site: https://<userid><at>stash.csiro.au/scm/~bow355/harman_r.git
 * 
 * Class CSelectRandom - used for random selection of K values for N, with or without replacement.
 */

enum replacement_type {WITHREPLACEMENT, WITHOUTREPLACEMENT} ;


#include <cstdlib>
#include <ctime>

// From: https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html
// Most of the supplied uniform generators return 32-bit integer values that are converted to doubles, 
// so they take at most 2^32 distinct values and long runs will return duplicated values
#if (_USE_RCPP==1)
#define _RANDFUNCTION_ (unif_rand() * (4294967295ul))
// N.B. #define _RANDSEED_FUNCTION will not do anything for the R RNG
#define _RANDSEED_FUNCTION ;
//#define _RANDFUNCTION_ rand()
#define _PRINTERROR Rcpp::Rcerr
#define _PRINTSTD Rcpp::Rcout
#define _USE_OMP 0
#else
#define _RANDFUNCTION_ rand()
#define _RANDSEED_FUNCTION srand(seedIN)
#define _PRINTERROR std::cerr
#define _PRINTSTD std::cout
#define _USE_OMP 1
#endif

// #define _SEEDCOUNTER std::chrono::system_clock::now().time_since_epoch().count()
// #define _SEEDCOUNTER GetCycleCount() 
#define _SEEDCOUNTER time(NULL)


template <class T>
class CSelectRandom 
{
	size_t  _CURRENTSEED ;
	bool    SEEDISSET ;
	void    SETSEED(size_t seedIN) ;

	public :

	CSelectRandom( size_t seedIN ) ;
	CSelectRandom( void) ;

	size_t GETSEED( void ) ;
	void FORCENEWSEED(size_t seedIN) ;

	std::vector<T>			        * Select(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN, replacement_type Replace, size_t rand_seed,  bool checkBool) ;
	std::vector<T>			        * SelectWITHREPLACEMENT(std::vector<T> * v_inputDataIN, size_t selectK_IN, size_t num_repeatsIN,  std::vector<size_t> * v_random64bitRangeIN) ;
	std::vector<T>			        * SelectWITHREPLACEMENTReturnRowAve(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN, std::vector<size_t> * v_random64bitRangeIN) ;

	std::vector<size_t>         *	ReturnVectOf64bitIntegersInRange_CSTDLIB(size_t howManyIN, size_t rand_seedIN, size_t startRangeIN, size_t endRangeIN) ;
	size_t					            GetRandom64Bits (unsigned int numbitswantedIN, unsigned int numbitsinrandIN ) ;

	std::vector<unsigned int>   *	ReturnVectOf32bitIntegersInRange_CSTDLIB(size_t howManyIN, size_t rand_seedIN, size_t startRangeIN, size_t endRangeIN) ;
	unsigned int 			          GetRandom32Bits (unsigned int numbitswantedIN, unsigned int numbitsinrandIN ) ;


  size_t 					            GetNumbitsRequired ( size_t inputData )  ;
  size_t					            GetNumRandomBitsRequired ( size_t maxBitsNeededIN, size_t randSizeIN )   ;
	//size_t					          GetRandom64Bits ( int numbitswanted ) ;


	std::vector<T> *		        SelectWITHOUTREPLACEMENT(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN) ;
	std::vector<T> *		        SelectWITHOUTREPLACEMENTReturnRowAve(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN) ;
	std::vector<T> *		        SelectWITHOUTREPLACEMENTReturnRowSumFast(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN) ;
	std::vector<T> * 		        SelectWITHOUTREPLACEMENTReturnRowSumFastTEST(size_t SizeOfInput, size_t select_KIN, size_t num_repeatsIN) ;

	void					         TestDistributionAverageRows( std::vector<T> * v_inputDataIN, T * permsDataIN, size_t selectK_IN, size_t num_repeatsIN  );

} ; // end of template <T> class CSelectRandom




template <class T>
 CSelectRandom<T>::CSelectRandom( void )
 {
	this->SEEDISSET = false ;
	this->_CURRENTSEED = 0 ;
 }

template <class T>
 CSelectRandom<T>::CSelectRandom( size_t seedIN )
 {
	if (seedIN == 0)
		seedIN = _SEEDCOUNTER ;

	this->SETSEED(seedIN) ;  // sets this->SEEDISSET == TRUE
	this->_CURRENTSEED = seedIN ;
 }



/*!  Select() - Selects data from an array at random specifying (WITHREPLACEMENT | WITHOUTREPLACEMENT)
	* \param v_inputDataIN   - Array of Y data (y is the response data)
	* \param selectK_IN			 - The number of values to select from the v_inputDataIN data (with or without replacement)
	* \param num_repeatsIN   - The number of permutations. The number of times to select 'num_select_cols' values from the original array 
	* \param Replace         - enum WITHREPLACEMENT, WITHOUTREPLACEMENT The selected values are to be repalced if \code true
	* \param rand_seed       - The random seed to use to initialise the random number generator
	* \param checkBool       - A flag to do a check of the algorithm. Requires num_select == size_y for a correct test. Also implements the timing of the function
	* \comment The CPU version has data in column major order as defined in original mpMap code
	*/
template <class T>
std::vector<T> * CSelectRandom<T>::Select(std::vector<T> * v_inputDataIN, size_t selectK_IN, size_t num_repeatsIN, replacement_type Replace, size_t rand_seed, bool checkBool)
{ 
		// size_t num_select  is the number of values to randomly select (number of cols in output matrix)
		// size_t nindividuals_cols = size_y ; size_y is the number of values in the Y vector
		size_t size_y = v_inputDataIN->size() ;
		if ((selectK_IN > size_y)&&(Replace == WITHOUTREPLACEMENT)) {
			_PRINTERROR << "Select() Error: Sellecting more values from a vector (without replacement) than what it contains" << std::endl ;
			return (NULL) ;
		}
    
		size_t howManyIN ;
		howManyIN = selectK_IN * num_repeatsIN ;

		std::vector<size_t> * v_randomArrayIndex = ReturnVectOf64bitIntegersInRange_CSTDLIB(howManyIN, rand_seed,  0, size_y-1 ) ;

		std::vector<T> * returnRandoms ;
		if ( Replace == WITHREPLACEMENT) 
		{
			returnRandoms = SelectWITHREPLACEMENT(v_inputDataIN,  selectK_IN,  num_repeatsIN, v_randomArrayIndex) ;
		}
		else if (Replace == WITHOUTREPLACEMENT)
		{
		//	returnRandoms = SelectWITHOUTREPLACEMENT(v_inputDataIN,   selectK_IN,  num_repeatsIN, v_randomArrayIndex ) ;
		}
    
  
	return (returnRandoms) ; // size will be [num_perms]x[num_select]
}  // end of sample()



/*!  SelectSelectWITHREPLACEMENT()
	* \param v_inputDataIN        - Vector of input values to select from
	* \param selectK_IN						- The number of values to select from the y data (with or without replacement)
	* \param num_repeatsIN				- The number of permutations. The number of times to select 'select_K' values from the original array 
	* \param v_random64bitRangeIN - vector of indexes into the inputData vector
	* \return result_perms				- an array of dimension [num_repeatsIN x selectK_IN] containg randomly selected elements of Y_Input[] values
	* \comment The CPU version has data in column major order as defined in original mpMap code
	*/
template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHREPLACEMENT(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN, std::vector<size_t> * v_random64bitRangeIN)
{

		// num_perms = the number of permutations to create - each one stored as a row of the output matrix
		T * ptr_inputData = v_inputDataIN->data() ;
				
			// This is the output array
		std::vector<T> * result_perms = (std::vector<T> *) new std::vector<T>(num_repeatsIN * selectK_IN) ;
		T * ptr_res_perms = (T*) result_perms->data() ;
		if (ptr_res_perms == NULL) {
			_PRINTERROR << "SelectWITHREPLACEMENT() Error: could not allocate results vector" << std::endl ;
			return (NULL) ;
		}

		size_t * ptr_randoms = v_random64bitRangeIN->data() ; // the array of random numbers between 0.0..1.0

		T newRand ;
		for (size_t h = 0 ; h < num_repeatsIN ; h++)
		{
			size_t currindex = h * selectK_IN ;
			for (size_t i = 0 ; i < selectK_IN; i++)
			{
				//newRand = distribution(rand_unif_generator) ;
				newRand = ptr_randoms[currindex + i] ; // get random value 				
				ptr_res_perms[currindex +i] = ptr_inputData[newRand]  ; // This part does the data selection
			}
		}
    
		return (result_perms) ; // size will be [num_repeatsIN]x[num_select]
}  // end of sample with replacement



 
/*!  SelectWITHREPLACEMENTReturnRowAve()
	* \param v_inputDataIN    - Array of data to be selected from (y is the response data)
	* \param selectK_IN				- The number of values to select from the y data (with or without replacement)
	* \param num_repeatsIN		- The number of permutations. The number of times to select 'select_K' values from the original array
	* \return result_perms		- an array of dimension [num_repeatsIN x selectK_IN] containg randomly selected elements of Y_Input[] values
	* \comment The CPU version has data in column major order as defined in original mpMap code
	*/
template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHREPLACEMENTReturnRowAve(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN, std::vector<size_t> * v_random64bitRangeIN)
{

		// num_perms = the number of permutations to create - each one stored as a row of the output matrix
	//	size_t size_vect  = v_inputDataIN->size() ;
		T * ptr_inputData = v_inputDataIN->data() ;
		std::vector<T> * result_perms = NULL ;

		if (selectK_IN > 0)
		{
				// This is the output array
				result_perms = (std::vector<T> *) new std::vector<T>(num_repeatsIN) ;
				T * ptr_res_perms = (T*) result_perms->data() ;
				if (ptr_res_perms == NULL) 
				{
					_PRINTERROR << "SelectWITHREPLACEMENT() Error: could not allocate results vector" << std::endl ;
					return (NULL) ;
				}

				size_t * ptr_sizet_randoms = v_random64bitRangeIN->data() ;

				size_t i_newRand ;
				T temp_ave ;
				for (size_t h = 0 ; h < num_repeatsIN ; h++)
				{
					//	nidivids = (T) size_vect - 1.0 ;  // this subtraction may cause issues if 
					temp_ave = 0.0 ; 
					for (size_t i = 0 ; i < selectK_IN; i++)
					{
						i_newRand = ptr_sizet_randoms[h * selectK_IN + i] ; // get random value 
						temp_ave += ptr_inputData[i_newRand]  ; // This part does the data selection
					}
					ptr_res_perms[h] = temp_ave / (double) selectK_IN ;
				}
		}
		else
		{
				result_perms = (std::vector<T> *) new std::vector<T>(0) ;
		}


    
		return (result_perms) ; // size will be [num_repeatsIN]x[num_select]
}  // end of sample with replacement



/*! \brief Function to randomly choose K values from a vector of data a desired number of times, so uses selection without replacement.
	* \param v_inputDataIN  - the data vector to select from
	* \param selectK_IN     - how many scalars to select from the vector
	* \param num_repeatsIN  - how many time to repeat the selection of K from the vector
	* \return A vector of size (num_repeatsIN * selectK_IN) that holds the selected values that were choosen from the inputData vector
	* \comment 
	*/
// Could think about using std::random_shuffle()
template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHOUTREPLACEMENT(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN)
{
		
		T * Y_Input = v_inputDataIN->data() ;
		std::vector<T> * ret_arrayOfSelections = NULL ; // the returned selected data - all K selctions for each repeat


		// size_t num_select  is the number of values to randomly select (number of cols in output matrix)
		size_t size_y ; 
		size_y = v_inputDataIN->size() ;

		if ((selectK_IN > size_y)) {
			_PRINTERROR << "Select() Error: Selecting more values from a vector (without replacement) than what it contains" << std::endl ;
			return (NULL) ;
		}


		if (selectK_IN > 0) 
		{
			
			// Copy input vector to muliplte std::list 
			//std::vector<std::list<T> > vector_of_y_as_list(num_repeatsIN) ;
			//for (size_t i = 0 ; i < num_repeatsIN; i++)
			//		vector_of_y_as_list[i](Y_Input, Y_Input+size_y);  // or could use: std::copy( v_inputDataIN->begin(), v_inputDataIN->end(), std::back_inserter( list ) );
			std::list<T> * inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;

			//std::copy( v_inputDataIN->begin(), v_inputDataIN->end(), inputData_as_list );

			// Create arrays of indexes that reduce in range for use as we select from an array
			// without replacement (i.e. selection from an array (list) that is shrinking after each selection)
			// This could cause some memory blowout on bigger num_repeats and selectK values.
			size_t nidivids ;
			size_t seed ;
			seed = 0 ; // if 0 then uses cueernt time as seed
			nidivids = size_y  ;
			std::vector<std::vector<size_t> *> t_arraysOfOffsets ;
			std::vector<size_t> * t_vec_randsinRange ;
			for (size_t i = 0 ; i < selectK_IN; i++)
			{
				t_vec_randsinRange = this->ReturnVectOf64bitIntegersInRange_CSTDLIB(num_repeatsIN, 0, 0, nidivids) ;
				t_arraysOfOffsets.push_back(t_vec_randsinRange) ;
				nidivids-- ;
			}
			
			// This is the output array 
			ret_arrayOfSelections= (std::vector<T> *) new std::vector<T>(num_repeatsIN * selectK_IN) ;
			if (ret_arrayOfSelections == NULL) {
				_PRINTERROR << "Select() Error: could not allocate results vector" << std::endl ;
				return (NULL) ;
			}
  
     size_t rand_offset ;
		 size_t currindex ;
			std::list<double>::const_iterator  it1 ;
			// This part does the data selection as goverend by the pre-computed offsets (arrayOfOffset) values.
			// Each selected element is removed from the list, before the next element is chosen.
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{

				currindex = h * selectK_IN ;
				for (size_t i = 0 ; i < selectK_IN; i++)
				{       
					it1 = inputData_as_list->begin() ;
					rand_offset = t_arraysOfOffsets[selectK_IN]->at(h) ;
					std::advance(it1,rand_offset) ; // this is the slow bit - advancing through a list is linear in the number of elements to traverse
					ret_arrayOfSelections->at(currindex +i)= *it1  ;  // save the randomly selected value
					it1 = inputData_as_list->erase(it1) ;  // remove the value that was selected (constant time)
				}

				delete inputData_as_list  ;
				inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			}

			delete inputData_as_list ;
			// cleanup array indexes
			for (size_t i = 0 ; i < selectK_IN; i++)
				delete t_arraysOfOffsets[i]  ;

		}
		else
			ret_arrayOfSelections = (std::vector<T> *) new std::vector<T>(0) ;
  

		return (ret_arrayOfSelections) ; 

    
	// _mm_free(arrayOfOffset) ;
	// size will be [num_perms]x[num_select]
}  // end of sample without replacement


template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHOUTREPLACEMENTReturnRowAve(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN)
{
		
		T * Y_Input = v_inputDataIN->data() ;
		std::vector<T> * ret_arrayOfSelections = NULL ; // the returned selected data - all averages - one for each repeat
		// size_t num_select  is the number of values to r andomly select (number of cols in output matrix)
		size_t size_y ; 
		size_y = v_inputDataIN->size() ;

		if ((selectK_IN > size_y)) {
			_PRINTERROR << "Select() Error: Selecting more values from a vector (without replacement) than what it contains" << std::endl ;
			return (NULL) ;
		}

		if (selectK_IN > 0) 
		{
			
			// Copy input vector to muliplte std::list 
			//std::vector<std::list<T> > vector_of_y_as_list(num_repeatsIN) ;
			//for (size_t i = 0 ; i < num_repeatsIN; i++)
			//		vector_of_y_as_list[i](Y_Input, Y_Input+size_y);  // or could use: std::copy( v_inputDataIN->begin(), v_inputDataIN->end(), std::back_inserter( list ) );
			// std::list<T> * inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			 std::vector<T> * inputData_as_list = (std::vector<T> *) new std::vector<T>(Y_Input, Y_Input+size_y) ;


			// Create arrays of indexes that reduce in range for use as we select from an array
			// without replacement (i.e. selection from an array (list) that is shrinking after each selection)
			// This could cause some memory blowout on bigger num_repeats and selectK values.
			size_t nidivids ;
			size_t seed ;
			seed = 0 ; // if 0 then uses cueernt time as seed
			nidivids = size_y  ;
			std::vector<std::vector<size_t> *> t_arraysOfOffsets ;
			std::vector<size_t> * t_vec_randsinRange ;
			for (size_t i = 0 ; i < selectK_IN; i++)
			{
				t_vec_randsinRange = this->ReturnVectOf64bitIntegersInRange_CSTDLIB(num_repeatsIN, 0, 0, nidivids) ;
				t_arraysOfOffsets.push_back(t_vec_randsinRange) ;
				nidivids-- ;
			}
			
			// This is the output array 
			ret_arrayOfSelections= (std::vector<T> *) new std::vector<T>(num_repeatsIN) ;
			if (ret_arrayOfSelections == NULL) {
				_PRINTERROR << "Select() Error: could not allocate results vector" << std::endl ;
				return (NULL) ;
			}
  
     size_t rand_offset ;
		 size_t currindex ;
		 T rowAve ;
			std::vector<double>::iterator  it1, it2 ;
			// This part does the data selection as goverend by the pre-computed offsets (arrayOfOffset) values.
			// Each selected element is removed from the list, before the next element is chosen.
	/*		This is *really really really* slow - uses a list of the input data.
			std::list<T>::iterator  it1, it2 ;
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{

				//currindex = h * selectK_IN ;
				rowAve = 0.0 ;
				for (size_t i = 0 ; i < selectK_IN; i++)
				{       
					it1 = inputData_as_list->begin() ;
					rand_offset = t_arraysOfOffsets[i]->at(h) ;
					std::advance(it1,rand_offset) ; // this is the slow bit - advancing through a list is linear in the number of elements to traverse
					rowAve += *it1  ;  // save the randomly selected value
					it1 = inputData_as_list->erase(it1) ;  // remove the value that was selected (constant time)
				}
				ret_arrayOfSelections->at(h) = rowAve / (T) selectK_IN ;

				delete inputData_as_list  ;
				inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			} */

			/* This works but is still slow - use std::vector of input data and std::vector<T>::iterators
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{

				//currindex = h * selectK_IN ;
				rowAve = 0.0 ;
				it2 = inputData_as_list->end() ;
				for (size_t i = 0 ; i < selectK_IN; i++)
				{       
					it1 = inputData_as_list->begin() ;					
					--it2 ;
					// std::advance(it2,-(i+1)) ;
					rand_offset = t_arraysOfOffsets[i]->at(h) ;
					std::advance(it1,rand_offset) ; // this is the slow bit - advancing through a list is linear in the number of elements to traverse
					rowAve += *it1  ;  // save the randomly selected value
					std::swap(*it1, *it2);  // swap the selected data with the data at the end of the list

				//	it1 = inputData_as_list->erase(it1) ;  // remove the value that was selected (constant time)
				}
				ret_arrayOfSelections->at(h) = rowAve / (T) selectK_IN ;

			//	delete inputData_as_list  ;
			//	inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			} */

			T lastval, selectval ;
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{
				//currindex = h * selectK_IN ;
				rowAve = 0.0 ;
				
				for (size_t i = 0 ; i < selectK_IN; i++)
				{       
					lastval = inputData_as_list->at(size_y-1-i)  ;
					rand_offset = t_arraysOfOffsets[i]->at(h) ;
					selectval = inputData_as_list->at(rand_offset) ; // this is the slow bit - advancing through a list is linear in the number of elements to traverse
					rowAve += selectval  ;  // save the randomly selected value
					inputData_as_list->at(rand_offset) = lastval ;  // swap the selected data with the data at the end of the list
					inputData_as_list->at(size_y-1-i)  = selectval ;
				//	it1 = inputData_as_list->erase(it1) ;  // remove the value that was selected (constant time)
				}
				ret_arrayOfSelections->at(h) = rowAve / (T) selectK_IN ;

			//	delete inputData_as_list  ;
			//	inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			}


			delete inputData_as_list ;
			// cleanup array indexes
			for (size_t i = 0 ; i < selectK_IN; i++)
				delete t_arraysOfOffsets[i]  ;
		}
		else
			ret_arrayOfSelections = (std::vector<T> *) new std::vector<T>(0) ;
  

		return (ret_arrayOfSelections) ; 
    
	// _mm_free(arrayOfOffset) ;
	// size will be [num_perms]x[num_select]
}  // end of sample without replacement



template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHOUTREPLACEMENTReturnRowSumFast(std::vector<T> * v_inputDataIN,  size_t selectK_IN, size_t num_repeatsIN)
{
		
		T * selectData_ptr = v_inputDataIN->data() ;
		std::vector<T> * ret_arrayOfSelections = NULL ; // the returned selected data - all averages - one for each repeat
		// size_t num_select  is the number of values to randomly select (number of cols in output matrix)
		size_t size_y ; 
		size_y = v_inputDataIN->size() ;

		if ((selectK_IN > size_y)) {
			_PRINTERROR << "Select() Error: Selecting more values from a vector (without replacement) than what it contains" << std::endl ;
			return (NULL) ;
		}

		if (selectK_IN > 0) 
		{
		  
#if (_PRINT==1)
		  _PRINTSTD << "CSelectRandom() obout to compute ReturnVectOf32bitIntegersInRange vector, size =  " << selectK_IN << std::endl  ;
#endif
			
			// Copy input vector to muliplte std::list 
			//std::vector<std::list<T> > vector_of_y_as_list(num_repeatsIN) ;
			//for (size_t i = 0 ; i < num_repeatsIN; i++)
			//		vector_of_y_as_list[i](Y_Input, Y_Input+size_y);  // or could use: std::copy( v_inputDataIN->begin(), v_inputDataIN->end(), std::back_inserter( list ) );
			// std::list<T> * inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;

			// Create arrays of indexes that reduce in range for use as we select from an array
			// without replacement (i.e. selection from an array (list) that is shrinking after each selection)
			// This could cause some memory blowout on bigger num_repeats and selectK values.
			size_t nidivids ;
			size_t seed ;
			seed = 0 ; // if 0 then uses cueernt time as seed
			nidivids = size_y  ;
			std::vector<std::vector<unsigned int> *> t_arraysOfOffsets ;
			std::vector<unsigned int> * t_vec_randsinRange ;
			for (size_t i = 0 ; i < selectK_IN; i++)
			{
				t_vec_randsinRange = this->ReturnVectOf32bitIntegersInRange_CSTDLIB(num_repeatsIN, seed, 0, nidivids) ;
				t_arraysOfOffsets.push_back(t_vec_randsinRange) ;
				nidivids-- ;
			}
#if (_PRINT==1)
			_PRINTSTD << "CSelectRandom() obtained ReturnVectOf32bitIntegersInRange vector " << std::endl  ;
#endif
			// This is the output array 
			ret_arrayOfSelections= (std::vector<T> *) new std::vector<T>(num_repeatsIN) ;
			if (ret_arrayOfSelections == NULL) {
				_PRINTERROR << "Select() Error: could not allocate results vector" << std::endl ;
				return (NULL) ;
			}
  
			unsigned int rand_offset ;
			unsigned int current_max  ;
			T rowSum, temp_max, temp_select ;
			//std::list<T>::const_iterator  it1 ;
			// This part does the data selection as governed by the pre-computed offsets (t_arraysOfOffsets) values.
			// Gamble-Durstenfeld-Fisher-Yates algorithm from  stackoverflow.com
			// Each selected element is placed at the back of the list, before the next element is chosen from a reduced size (by 1) array.
			// 
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{
				//currindex = h * selectK_IN ;
				rowSum = 0.0 ;
				current_max = size_y - 1;
				for (size_t i = 0 ; i < selectK_IN; i++)
				{       					
					rand_offset = t_arraysOfOffsets[i]->at(h) ;  // the offsets are computed to be in the correct range of size of the (potentially reduced size) array
					temp_select = selectData_ptr[rand_offset] ;
					rowSum += temp_select ;
					temp_max = selectData_ptr[current_max] ;  // make a copy of current value from the back of the array
					selectData_ptr[current_max] = temp_select ;  // place selection at back of array
					selectData_ptr[rand_offset] = temp_max ;  // place original last value at the random position
					--current_max ;  // reduce current max so that it does not include the already selected element (which was placed at the  back of the array)
				}
				ret_arrayOfSelections->at(h) = rowSum  ;

			//	delete inputData_as_list  ;
			//	inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			}

		//	delete inputData_as_list ;
			// cleanup array indexes
			for (size_t i = 0 ; i < selectK_IN; i++)
				delete t_arraysOfOffsets[i]  ;
		}
		else
			ret_arrayOfSelections = (std::vector<T> *) new std::vector<T>(0) ;
  

		return (ret_arrayOfSelections) ; 
    
	// _mm_free(arrayOfOffset) ;
	// size will be [num_perms]x[num_select]
}  // end of sample without replacement



template <class T>
std::vector<T> * CSelectRandom<T>::SelectWITHOUTREPLACEMENTReturnRowSumFastTEST(size_t SizeOfInput, size_t select_KIN, size_t num_repeatsIN)
{
		std::vector<T> * v_inputData ;
		v_inputData = (std::vector<T> *) new std::vector<T>(SizeOfInput) ;
		T * selectData_ptr = v_inputData->data() ;
		T vectorSumRes ; // use this to check accuracy of results
		for (size_t t1 = 0 ; t1 < SizeOfInput; t1++)
		{
			selectData_ptr[t1] = (T)( t1 + 1.0 );
			vectorSumRes +=  t1 + 1.0 ;
		}

		std::vector<T> * ret_arrayOfSelections = NULL ; // the returned selected data - all averages - one for each repeat
		// size_t num_select  is the number of values to randomly select (number of cols in output matrix)
		size_t size_y ;
		size_y = v_inputData->size() ;
		unsigned int  fullRange  = size_y ;

		if ((select_KIN > size_y)) {
			_PRINTERROR << "Select() Error: Selecting more values from a vector (without replacement) than what it contains" << std::endl ;
			return (NULL) ;
		}

		if (select_KIN > 0)
		{

			// Copy input vector to muliplte std::list
			//std::vector<std::list<T> > vector_of_y_as_list(num_repeatsIN) ;
			//for (size_t i = 0 ; i < num_repeatsIN; i++)
			//		vector_of_y_as_list[i](Y_Input, Y_Input+size_y);  // or could use: std::copy( v_inputDataIN->begin(), v_inputDataIN->end(), std::back_inserter( list ) );
			// std::list<T> * inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;

			// Create arrays of indexes that reduce in range for use as we select from an array
			// without replacement (i.e. selection from an array (list) that is shrinking after each selection)
			// This could cause some memory blowout on bigger num_repeats and selectK values.
			size_t nidivids ;
			size_t seed ;
			seed = 0 ; // if 0 then uses cueernt time as seed
			nidivids = size_y  ;
			std::vector<std::vector<unsigned int> *> t_arraysOfOffsets ;
			std::vector<unsigned int> * t_vec_randsinRange ;
			for (size_t i = 0 ; i < select_KIN; i++)
			{
				t_vec_randsinRange = this->ReturnVectOf32bitIntegersInRange_CSTDLIB(num_repeatsIN, seed, 0, nidivids) ;
				t_arraysOfOffsets.push_back(t_vec_randsinRange) ;
				nidivids-- ;
			}

			// This is the output array
			ret_arrayOfSelections= (std::vector<T> *) new std::vector<T>(num_repeatsIN) ;
			if (ret_arrayOfSelections == NULL) {
				_PRINTERROR << "Select() Error: could not allocate results vector" << std::endl ;
				return (NULL) ;
			}

			struct timespec ts1, ts2 ;
			double time1, time2 ;

#if (TIMING)
			clock_gettime(CLOCK_MONOTONIC, &ts1) ;
#endif
			unsigned int rand_offset ;
			unsigned int current_max  ;
			T rowSum, temp_max, temp_select ;
			//std::list<T>::const_iterator  it1 ;
			// This part does the data selection as governed by the pre-computed offsets (t_arraysOfOffsets) values.
			// Gamble-Durstenfeld-Fisher-Yates algorithm from  stackoverflow.com
			// Each selected element is placed at the back of the list, before the next element is chosen from a reduced size (by 1) array.
			//
			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{
				//currindex = h * selectK_IN ;
				rowSum = 0.0 ;
				current_max = size_y - 1;
				for (size_t i = 0 ; i < select_KIN; i++)
				{
					rand_offset = t_arraysOfOffsets[i]->at(h) ;  // the offsets are computed to be in the correct range of size of the (potentially reduced size) array
					temp_select = selectData_ptr[rand_offset] ;
					rowSum += temp_select ;
					temp_max = selectData_ptr[current_max] ;  // make a copy of current value from the back of the array
					selectData_ptr[current_max] = temp_select ;  // place selection at back of array
					selectData_ptr[rand_offset] = temp_max ;  // place original last value at the random position
					--current_max ;  // reduce current max so that it does not include the already selected element (which was placed at the  back of the array)
				}
				ret_arrayOfSelections->at(h) = rowSum  ;

			//	delete inputData_as_list  ;
			//	inputData_as_list = (std::list<T> *) new std::list<T>(Y_Input, Y_Input+size_y) ;
			}


#if (TIMING)

			clock_gettime(CLOCK_MONOTONIC, &ts2) ;
			// get timing results
			time1 = (double) ts1.tv_sec * 1000000000.0 + (double) + (double) ts1.tv_nsec ;
			time2 = (double) ts2.tv_sec * 1000000000.0 + (double) + (double) ts2.tv_nsec ;
			_PRINTSTD << "Time for CPU selectNK: " << ( time2 - time1) / 1000000000.0 << std::endl ;

#endif


			_PRINTSTD << "CSelectRandom<T>::SelectWITHOUTREPLACEMENTReturnRowSumFastTEST Test results:" << std::endl ;
			_PRINTSTD << "Select: " << select_KIN << " From  " << fullRange << std::endl ;
			_PRINTSTD << "With input data range from 1.. " << SizeOfInput << std::endl ;
			size_t failedCount = 0 ;
			if (select_KIN == fullRange)
			{
				for (size_t t1 = 0 ; t1 < num_repeatsIN; t1++)
				{
						if (vectorSumRes != ret_arrayOfSelections->at(t1))
						{
							if (failedCount < 32)
							{
								_PRINTSTD << "Error at:\t\t" << t1 << "\tResult =\t" << ret_arrayOfSelections->at(t1) << "\t Expected:\t" << vectorSumRes << std::endl  ;
								failedCount++ ;
							}
							else
								failedCount++ ;
						}
				}
				// compute average:
				T aveTest = 0.0;
				for (size_t t1 = 0 ; t1 < num_repeatsIN; t1++)
					aveTest += ret_arrayOfSelections->at(t1) ;
				aveTest /= num_repeatsIN ;

				if (failedCount > 0)
					_PRINTSTD << "Failed test! Count of fails :\t" <<  failedCount << " \t from "<< num_repeatsIN <<  std::endl  ;
				else
					_PRINTSTD << "Passed test - All sums equaled expected result."  << std::endl  ;

				_PRINTSTD << "Permuted average: " << aveTest << " \t % expected:\t" << aveTest / vectorSumRes * 100<< "%" << std::endl ;
			}
			else
			{
				// compute average:
				T aveTest = 0.0;
				for (size_t t1 = 0 ; t1 < num_repeatsIN; t1++)
					aveTest += ret_arrayOfSelections->at(t1) ;
				aveTest /= num_repeatsIN ;

				_PRINTSTD << " Expected average:  " << (vectorSumRes / (T) fullRange)  << "\tPermuted average: " << aveTest/(T)select_KIN  << std::endl ;
			}
		//	delete inputData_as_list ;
			// cleanup array indexes
			for (size_t i = 0 ; i < select_KIN; i++)
				delete t_arraysOfOffsets[i]  ;
		}
		else
			ret_arrayOfSelections = (std::vector<T> *) new std::vector<T>(0) ;


		return (ret_arrayOfSelections) ;

	// _mm_free(arrayOfOffset) ;
	// size will be [num_perms]x[num_select]
}  // end of sample without replacement fastTEST


/************************************************************************
 * The stdlib rand() function will always be available:
 ************************************************************************/
// Use stdlib rand() function
template <class T>
size_t CSelectRandom<T>::GetRandom64Bits (unsigned int numbitswantedIN, unsigned int numbitsinrandIN )
{
#if (_USE_RCPP==1)
  //  GetRNGstate();
#if (_PRINT==2)
  static unsigned int got_here = 0 ;
  if (got_here == 0) {
    Rcpp::Rcout << "Using 64bit random function" << std::endl ; 
    got_here = 1 ;
  }
#endif
#endif
		size_t randData1, randData2 ;
		size_t gotbits = 0 ;
		randData1 = _RANDFUNCTION_ ;
		gotbits = numbitsinrandIN ;
		while (gotbits < numbitswantedIN)
		{		
			randData1 = randData1 << numbitsinrandIN ;
			gotbits += numbitsinrandIN ;
			randData2  = _RANDFUNCTION_ ;
			randData1 = randData1 | randData2 ;
		}
		
#if (_USE_RCPP==1)
//		PutRNGstate();
#endif
		
		return randData1 ;
}

template <class T>
unsigned int CSelectRandom<T>::GetRandom32Bits (unsigned int numbitswantedIN, unsigned int numbitsinrandIN )
{

		unsigned int randData1, randData2 ;
		unsigned int gotbits = 0 ;

		randData1 = _RANDFUNCTION_ ;
		gotbits = numbitsinrandIN ;
		while (gotbits < numbitswantedIN)
		{
			randData1 = randData1 << numbitsinrandIN ;
			gotbits += numbitsinrandIN ;
			randData2  = _RANDFUNCTION_ ;
			randData1 = randData1 | randData2 ;
		}
	
		return randData1 ;
}


// Use stdlib rand() function
template <class T>
size_t CSelectRandom<T>::GetNumbitsRequired ( size_t inputData ) 
{

    size_t numbitsinrand ;
		numbitsinrand = 0 ;
		while (inputData >= 1)
		{
			numbitsinrand++ ;
			inputData = inputData >> 1 ;
		}
#if (_PRINT==2)
		static unsigned int got_here = 0 ;
		if (got_here == 0) {
		  Rcpp::Rcout << "GetNumbitsRequired: " << numbitsinrand << std::endl ; 
		  got_here = 1 ;
		}
#endif
		return numbitsinrand ;
}



template <class T>
size_t CSelectRandom<T>::GetNumRandomBitsRequired ( size_t maxBitsNeededIN, size_t randSizeIN ) 
{
    assert(randSizeIN > 0) ;
    
    size_t retNumRands ;
  	retNumRands = 1 ;
	
    while ((retNumRands * randSizeIN) < maxBitsNeededIN)
		{
			retNumRands++ ;
		}
		
		return retNumRands  ;
}



template <class T>
void CSelectRandom<T>::SETSEED(size_t seedIN)
{
	if (seedIN == 0)
		seedIN = _SEEDCOUNTER ;

	_RANDSEED_FUNCTION ;
	//srand (seedIN);
	// std::default_random_engine rand_unif_generator(rand_seedIN)
	this->SEEDISSET = true ;


}

template <class T>
void CSelectRandom<T>::FORCENEWSEED(size_t seedIN)
{

  _RANDSEED_FUNCTION ;
	// srand (seedIN);
	// std::default_random_engine rand_unif_generator(rand_seedIN)
	this->SEEDISSET = true ;


}

template <class T>
size_t CSelectRandom<T>::GETSEED( void )
{
	// srand (seed);
	// std::default_random_engine rand_unif_generator(rand_seedIN)
	if (this->SEEDISSET == true)
		return this->_CURRENTSEED ;
	else
	{
		_PRINTERROR << "Error: CSelectRandom<T>::GETSEED(): seed has not been set" << std::endl ;
		return this->_CURRENTSEED ;
	}
}



template <class T>
std::vector<unsigned int> * CSelectRandom<T>::ReturnVectOf32bitIntegersInRange_CSTDLIB(size_t howManyIN, size_t rand_seedIN, size_t startRangeIN, size_t endRangeIN)
{

		if (this->SEEDISSET == false)
			this->SETSEED(rand_seedIN) ;


		std::vector<unsigned int> * v_randomDiscreteRangeData = NULL ;
		long long int fullRange  ;
		long long int cutoff ;
		fullRange = endRangeIN - startRangeIN ;

		if (fullRange > 0)
		{

		  long long int maxBitsNeeded ;
		  long long int randSizeBits ;
		  long long int numRandsNeeded ;
#if (_USE_RCPP==1)
			randSizeBits = this->GetNumbitsRequired(4294967295ul) ;  // we assume that the R docs here https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html are correct
#else
			randSizeBits = this->GetNumbitsRequired(RAND_MAX) ;	
#endif	
				
				maxBitsNeeded  = this->GetNumbitsRequired(fullRange) ;// get the number of bits that would represent the largest number needed
				numRandsNeeded = this->GetNumRandomBitsRequired(maxBitsNeeded, randSizeBits) ;

				unsigned int maxForNumbitsGiven ;
				long long int maxForNumbitsGiven_ll ;
				maxForNumbitsGiven_ll = 0 ;

				maxForNumbitsGiven_ll = ((long long int)1 <<   (unsigned int) (numRandsNeeded * randSizeBits ))  ;
				maxForNumbitsGiven_ll = maxForNumbitsGiven_ll - 1 ;
				maxForNumbitsGiven    = maxForNumbitsGiven_ll ;

				cutoff = (maxForNumbitsGiven / fullRange) * fullRange ;  // this is so we do not have a bias when setting the new range to the random values returned

				v_randomDiscreteRangeData = (std::vector<unsigned int> * ) new  std::vector<unsigned int>(howManyIN) ;
				unsigned int * ptr_dataOUT = v_randomDiscreteRangeData->data() ;

#if (_USE_RCPP==1)
#if (_PRINT==2)
				Rcpp::Rcout << "1:" << numRandsNeeded << "  maxBitsNeeded:" << maxBitsNeeded << ", maxForNumbitsGiven:" << maxForNumbitsGiven << ", cutoff:"<< cutoff << ", fullRange:" << fullRange << ", randSizeBits:" << randSizeBits << std::endl ;
#endif
#endif
				unsigned int randData1  ;
			//	size_t test_lli ; 
			long long int reject ;
				//size_t reject ;
				long long int currval ;
				reject = 0 ;
				currval = 0 ;
				//	for (size_t t1 = 0 ; t1 < v_random64bitIN->size() ; t1++)
#if (_USE_OMP==1)
#pragma omp critical 
#endif
{

				while (currval < howManyIN)
				{
					randData1 = this->GetRandom32Bits(maxBitsNeeded, randSizeBits) ;

					if (randData1 > cutoff)
						reject++ ;
					else
					{
						ptr_dataOUT[currval++] = startRangeIN + (randData1 % fullRange) ;
					}
				}
	
}
				v_randomDiscreteRangeData->resize(currval) ;
				return (v_randomDiscreteRangeData) ;
		}
		else if (fullRange == 0 )
		{
			return NULL ;
		}
		else // the range of values to select is 0, so no values can be selected from that set.
		{
			_PRINTERROR << "Error: CSelectRandom<T>::ReturnVectOf32bitIntegersInRange_CSTDLIB(): fullRange should be >= 0" << std::endl ;
			return NULL ;
		}

}


template <class T>
std::vector<size_t> * CSelectRandom<T>::ReturnVectOf64bitIntegersInRange_CSTDLIB(size_t howManyIN, size_t rand_seedIN, size_t startRangeIN, size_t endRangeIN) 
{

		if (this->SEEDISSET == false)
			this->SETSEED(rand_seedIN) ;

        
		std::vector<size_t> * v_randomDiscreteRangeData = NULL ;
		size_t fullRange  ;
		size_t cutoff ;   
		fullRange = endRangeIN - startRangeIN ;

		if (fullRange > 0)
		{

				size_t maxBitsNeeded ;
				size_t randSizeBits ;
				size_t numRandsNeeded ;
    
#if (_USE_RCPP==1)
    randSizeBits = this->GetNumbitsRequired(4294967295ul) ;  // we assume that the R docs here https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html are correct
#else
    randSizeBits = this->GetNumbitsRequired(RAND_MAX) ;	
#endif 
				maxBitsNeeded = this->GetNumbitsRequired(fullRange) ;// get the number of bits that would represent the largest number needed
				numRandsNeeded = this->GetNumRandomBitsRequired(maxBitsNeeded, randSizeBits) ;
    
				size_t maxForNumbitsGiven ;
				maxForNumbitsGiven = 0 ;

				maxForNumbitsGiven = ((size_t)1 <<   (size_t) (numRandsNeeded * randSizeBits ))  ;
				maxForNumbitsGiven -= 1 ;
				/*
			 	size_t temp_st = 1 ;
				for (size_t t1 = 0 ; t1 < (numRandsNeeded * randSizeBits ) ; t1++)
				{
					temp_st = temp_st << 1 ;
				}
				maxForNumbitsGiven = temp_st -1 ; */
				// std::numeric_limits<size_t>::max()  
				cutoff = (maxForNumbitsGiven / fullRange) * fullRange ;  // this is so we do not have a bias when setting the new range to the random values returned

				v_randomDiscreteRangeData = (std::vector<size_t> * ) new  std::vector<size_t>(howManyIN) ;
				size_t * ptr_dataOUT = v_randomDiscreteRangeData->data() ;
				
		
				size_t randData1, randData2  ;

				size_t reject ;
				size_t currval ;
				reject = 0 ;
				currval = 0 ; 
				//	for (size_t t1 = 0 ; t1 < v_random64bitIN->size() ; t1++)
#if (_USE_RCPP==1)
				GetRNGstate();
#endif
				while (currval < howManyIN)
				{
					randData1 = this->GetRandom64Bits(maxBitsNeeded, randSizeBits) ;

					if (randData1 > cutoff)
						reject++ ;
					else
					{
						ptr_dataOUT[currval++] = startRangeIN + (randData1 % fullRange) ;
					}
				}
#if (_USE_RCPP==1)
				PutRNGstate();
#endif				
  
				v_randomDiscreteRangeData->resize(currval) ;
				return (v_randomDiscreteRangeData) ;
		}
		else if (fullRange == 0 )
		{
			return NULL ;
		}
		else // the range of values to select is 0, so no values can be selected from that set.
		{
			_PRINTERROR << "Error: CSelectRandom<T>::ReturnVectOf64bitIntegersInRange_CSTDLIB(): fullRange should be >= 0" << std::endl ;
			return NULL ;
		}

}
//#endif




/*!  
	* \param inputDataIN     - Array of input data that was selected from
	* \param permsDataIN      - The size of the Y_Input array (in elements of type T) 
	* \param num_repeatsIN		- The number of permutations. The number of times to select 'select_K' values from the original array
	* \return result_perms		- an array of dimension [num_repeatsIN x selectK_IN] containg randomly selected elements of Y_Input[] values
	* \comment The CPU version has data in column major order as defined in original mpMap code
*/
template <class T>
void CSelectRandom<T>::TestDistributionAverageRows( std::vector<T> * v_inputDataIN, T * permsDataIN, size_t selectK_IN,  size_t num_repeatsIN  ) 
{
		size_t size_vect = v_inputDataIN->size() ;
	// test the output data:
	//	if (checkBool == true)
		{
			_PRINTSTD << "The average of  the resampled row values should approach the average of the  original data at high select_K" << std::endl ;
			size_t numMore = 0; 
			size_t numLess = 0;
			double inputAve = 0.0 ;
			T * ptr_inputData = v_inputDataIN->data() ;
			for (size_t i = 0 ; i < size_vect; i++) 
				inputAve += ptr_inputData[i] ;
			inputAve = inputAve / (double) size_vect  ;

			for (size_t h = 0 ; h < num_repeatsIN ; h++)
			{
				size_t currindex = h * selectK_IN ;
				double rowave = 0.0 ;
				for (size_t i = 0 ; i < selectK_IN; i++)
				{     
					rowave += permsDataIN[currindex + i]   ;              
				} 
				rowave = rowave / (double) selectK_IN ;
				// the average of  the resampled row values should be approach to the average of the  original data at high selectK_IN
				// _PRINTSTD <<  " row: "<< h << " sum: "<< rowave << " expected: "<< firstRowAve << " diff: " << firstRowAve - rowave <<  std::endl ;
				if (inputAve - rowave > 0.0) 
					numMore++ ;
				else 
					numLess++ ;
				}
				_PRINTSTD <<  " numMore: "<< numMore << " numLess: "<< numLess << " These should be approximately equal!" <<  std::endl ;
			}
}

 // end of implementation of template <T> class CSelectRandom
