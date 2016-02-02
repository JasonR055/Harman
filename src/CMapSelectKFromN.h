/**
 * \file Harman/src/CMapSelectKFromN.h 
 * \date 12/11/2015
 * \author   Josh Bowden, CSIRO
 * \details Harman eResearch Collaboration project ERRFP-263 (https://jira.csiro.au/browse/ERRFP-263)
 * Stash site: https://<userid><at>stash.csiro.au/scm/~bow355/harman_r.git
 * 
 * Class CFactorials:  Computes the factorial of positive input integers and stores them as an array for later access
 * Class CMapSelectKFromN: Creates a map of the scans (progressive sum of a vector) of all recursively defined binomial 
 * sequences relating to choosing K values from the set of N. The functions GetScan_X() allow the caller to return a set of K 
 * indexes into an array of N values. The set of K values is dependent on the input value (which ranges from 0-nchhosek(N,K)) 
 * and they correspond to a distict order of arranging the possible ways to return K values of a set of N. 
 * e.g. [0,1,2,3,4,5,6,7] N = 8, choose K = 4
 * 0,1,2,3
 * 0,1,2,4
 * 0,1,2,5
 * 0,1,2,6
 * 0,1,2,7
 * 0,1,3,4
 * 0,1,3,5
 * 0,1,3,6
 * 0,1,3,7
 * 0,1,4,5
 * 0,1,4,6
 * 0,1,4,7
 * ... etc.
 */

#include <string>

#if (_USE_RCPP==1)
#define _PRINTERROR Rcpp::Rcerr
#define _PRINTSTD Rcpp::Rcout
#else
#define _PRINTERROR std::cerr
#define _PRINTSTD std::cout
#endif

/*! \brief Class that computes the factorial of positive input integers and stores them as an array for later access
 *  \reference 
 * 
 */
//template <class T> class CFactorials 

template <class T> class CFactorials
{
private:
  std::vector<T> v_factorialData ;
 
public:
  // Constructor creates an array of factorial results for numbers 1..maxN
  CFactorials<T>(size_t maxN)
  { 
    this->v_factorialData.reserve(32) ; // reserve space for all calculated values
    
    // set first value
    this->v_factorialData.push_back((T)1) ;
    
    // calculate the factorial values iteratively
    for(size_t i = 1; i <= maxN; i++)
    {
      // this->v_factorialData[i] = this->v_factorialData[i-1] * i ;
      this->v_factorialData.push_back(this->v_factorialData[i-1] * i) ;
    }
          
  }
  
  /*! Returns the factorial of the input value */
  T nfact(size_t nIN)
  {  
      // make sure result has already been calculated
      if (nIN < this->v_factorialData.size() )
      {         
        return this->v_factorialData[nIN];
      }
      else
      {
        UpadateFactorialData(nIN) ;
        return this->v_factorialData[nIN];
      }
  }
  
	/*!
   * Extends the v_factorialData std::vector<> so that it contains the factorial values up to the input nIN 
   */
  void UpadateFactorialData(size_t nIN)
  {
    if (nIN > 20) {
      if (typeid(T) == typeid(size_t))
      {
        std::cout << std::flush ;
        _PRINTERROR << "StoreFactorial() Error: maximum factorial in 64 bit integer is 20!. Truncating input." << std::endl ;
        nIN = 20 ;
      }
    }
    
    size_t old_size = this->v_factorialData.size() ;
    this->v_factorialData.reserve(nIN) ; // reserve space for all calculated values
    
    // calculate the factorial values iteratively
    for(size_t i = old_size; i <= nIN; i++)
        this->v_factorialData.push_back(this->v_factorialData[i-1] * i) ;

  }
  
  /*!
   * Prints values in the std::vector<> containing the factorial values to the screen 
   */
  void Print( void )
  {
    for(size_t i = 0 ; i < this->v_factorialData.size(); i++)
          std::cout <<  this->v_factorialData[i] << " "  ;
    std::cout << std::endl ;
  }
  
  
  /*! An iterative method to calculate the factorial of a number */
  T iter_factorial(size_t nIN)
  {
      T ret = (T) 1;
      for(size_t i = 1; i <= nIN; ++i)
          ret *= (T) i;
      return ret;
  }
  
  /*! Calculate a partial factorial : i.e. (8,4) = 8 x 7 x 6 x 5 
  * useful for the 
  */
  T iter_factorial_numerator(size_t nIN, size_t downtoIN)
  {
      T ret = (T) 1;
      for(size_t i = nIN; i > downtoIN ; i--)
          ret *= (T) i;
      return ret;
  }
  
  
  
  /*! Calculates the binomial coeficient of (N,K)
	*
	* @description
  * \code{nchoosek} Calculates the binomial coeficient of (N,K) - total number of possiblities of from set of N values, select K of them
	*  It should keep intermediate results within the numeric accuaracy of 64 bit computers until n-k is large
  *  This version chooses the least number of divisions.
  */ 
  T nchoosek(size_t nIN, size_t kIN)
  { 
    
    T ret = 1.0 ;
    T bigger  ;
    T smaller ;
    T nmink = (double)  nIN -  (double) kIN ;
    

		// makes use of the fact that (N,K) == (N,N-K)
		// e.g. (20,3) == (20,(20-3))
    if (nmink <  kIN) {
      bigger  = (double) kIN ;
      smaller = nmink ; 
    } else  {  // (nmink >=  kIN)
      bigger  = nmink ;
      smaller= (double)  kIN ;
    }    
    
/*    fpu_control_t _oldcw, _cw;
    _FPU_GETCW(_oldcw);
    _cw = (_oldcw & ~_FPU_DOUBLE & ~_FPU_SINGLE) | _FPU_EXTENDED;
    _FPU_SETCW(_cw);
  */
      //seq_nk <- seq(N,bigger+1, -1 )
      for(size_t i = nIN ; i > bigger; i--)
      {
        ret =  ret * ( (T) i / (T)  smaller) ;                
        smaller =  smaller - 1 ;
      }
      
 //     _FPU_SETCW(_oldcw);
      return ((T) ret);
  }

  
  /*
	*
  * It should keep intermediate results within the numeric accuaracy of 64 bit computers until n-k is  large
  */ 
 /* T nchoosek(size_t nIN, size_t kIN)
  {   
      assert(nIN >= kIN) ;
      assert(nIN >= 0) ;
      assert(kIN >= 0) ;
      
      double ret = (double) 1.0 ;
      size_t nmink = nIN - kIN ;
      
#if !defined (_MSC_VER)
      Rcpp::Rcout.precision(25);
      Rcpp::Rcout.setf( std::ios::fixed, std:: ios::floatfield );
#endif

     if (nIN > kIN)
     {
         for(size_t i = nIN ; i > kIN ; i--)
         {     
            if (nmink > 0)
              ret *= (double) i / (double) (nmink) ;
             
    #if defined (_MSC_VER)
              std::cout << " ret :\ti:"  << i << "\t nmink:" << nmink << "  \t"  <<  ret << std::endl ;  std::cout << std::flush ;
    #else 
              Rcpp::Rcout << " ret :\ti:"  << i << "\t nmink:" << nmink << "  \t"  <<  ret << std::endl ;  std::cout << std::flush ;
    #endif                 
             
              nmink-- ;
              if (nmink < 1) 
                nmink = 1 ;
                       
         }
     }
     // else if (nIN == kIN) return (1)
               
      return round(ret);
  }
 */

/*
   * \returns n! / (n-k)! x k!  = the number of combinations of selecting k values from a set of n values, except if n or k are ==0, when it returns 1.
   */
  /*
  T nchoosek(size_t nIN, size_t kIN) 
  {
    assert(nIN >= kIN) ;
    assert(nIN >= 0) ;
    assert(kIN >= 0) ;
    
    T nf =  this->nfact(nIN) ;
    T kf =  this->nfact(kIN) ;
    T nminuskf =  this->nfact(nIN-kIN) ;
    T res = nf / (nminuskf * kf) ;
    
    if (kIN == 0) res = (T) 1 ;  
    return (res) ;
  }
  */
/*  T nchoosek(size_t nIN, size_t kIN) 
  {
    assert(nIN >= kIN) ;
    assert(nIN >= 0) ;
    assert(kIN >= 0) ;
    
    T numerator = this->iter_factorial_numerator(nIN, kIN) ;  
	T nminuskf =  this->nfact(nIN-kIN) ; // hopefully this remains in our full precision range - It will not if the difference (n-k) is > 18
   // std::cout << " (n,k) = " << nIN<<  "," << kIN << ",\tnminuskf: " << nminuskf << ",\tnumerator: " << numerator  <<  std::endl ;   
    T res =  numerator / nminuskf ;
    
    if (kIN == 0) res = (T) 1 ;  
    return (res) ;
  }*/
  
   void TestFactorialClass() 
  {
    
	  std::cout.precision(15);
	  std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed    
    std::cout << "Running tests in CFactorials::TestFactorialClass()  " << std::endl ;
    
    std::cout  << " nchoosek(8,4) :\t"  	<<  this->nchoosek(8,4) << std::endl ;
    double ret = this->nchoosek(300,20) ;
    //std::cout << " nchoosek(300,20) :\t"  <<  ret   << " \tfrom R: 7.500434e+30" << std::endl ;
    std::cout << " got here ok:\t" << std::endl ;
    std::cout  << ret << std::endl ;
    std::cout << " got here ok:\t" << std::endl ;
    /*
    std::cout << " nchoosek(6,2) :\t"		<<  this->nchoosek(6,2) << std::endl ;
    std::cout << " nchoosek(8,4) :\t"		<<  this->nchoosek(8,4) << std::endl ;
    std::cout << " nchoosek(8,2) :\t"		<<  this->nchoosek(8,2) << std::endl ;
    std::cout << " nchoosek(4,2) :\t"		<<  this->nchoosek(4,2) << std::endl ;
    std::cout << " nchoosek(6,0) :\t"		<<  this->nchoosek(6,0) << std::endl ;
    std::cout << " nchoosek(0,0) :\t"		<<  this->nchoosek(0,0) << std::endl ;
   // std::cout << " nchoosek(0,4) :"	<<  this->nchoosek(0,4) << std::endl ;
    
    std::cout << " nchoosek(30,4) :\t"		<<  this->nchoosek(30,4)   << " \tfrom R: 27405" << std::endl ;
    std::cout << " nchoosek(20,10) :\t"		<<  this->nchoosek(20,10)  << " \tfrom R: 184756" << std::endl ;
    std::cout << " nchoosek(100,90) :\t"	<<  this->nchoosek(100,90) << " \tfrom R: 1.731031e+13" << std::endl ;
    std::cout << " nchoosek(50,10) :\t"		<<  this->nchoosek(50,10)  << " \tfrom R: 10272278170" << std::endl ;
	  std::cout << " nchoosek(26,7) :\t"		<<  this->nchoosek(26,7)   << " \tfrom R: 65800" << std::endl ;
    std::cout << " nchoosek(300,1) :\t"		<<  this->nchoosek(300,1)   << " \tfrom R: 300" << std::endl ;
	  std::cout << " nchoosek(300,20) :\t"	<<  this->nchoosek(300,20)   << " \tfrom R: 7.500434e+30" << std::endl ;

    std::cout << " i= " << 18 << " "	<<	this->nfact(18) << " from R:  6.402374e+15"	<< std::endl ; ;
    std::cout << " i= " << 20 << " "	<<	this->nfact(20) << " from R:  2.432902e+18"	<< std::endl ; ;
    std::cout << " i= " << 25 << " "	<<	this->nfact(25) << " from R:  1.551121e+25"	<< std::endl ; ;
    for (size_t i = 0 ; i < 24 ; i++)
      std::cout << " i= " << i << " " << this->nfact(i) << " " << std::endl ; ;
     */   
    std::cout << std::endl ;    
  }
  
 } ; // end of class CFactorials 




typedef std::pair<size_t, size_t> key_NK_pair ;  // requires #include <utility>
  
/*! \brief Creates a map of the scans (progressive sum of a vector) of all recursively defined binomial sequences relating to choosing K values from the set of N.
 * 
 *  \details Creates all binomial sequences for input numbers (N,K), N being the number in a set and K being the size of the sample to take form the set.
 *    Each sequence is the scan of the binomial coefficients (N-1,K-1),(N-2,K-1),(N-3,K-1)...(N-X,K-1) where X = N-K-1
 *    e.g. (8,4)  =  (7,3),(6,3),(5,3),(4,3),(3,3) where X = 8-4-1 = 3
 *                =   35  , 20,   10,   4,    1
 *          scan  =   35,   55,   65,  69,  70
 *    Results in a map of each sequence scan from N..1 and K..1
 * 
 */
class CMapSelectKFromN
{
private:

  std::map< key_NK_pair, std::vector<size_t >* > map_nchoosek ;  /*! This holds all the sequences related to the input (N,K) = (N-1,K-1),(N-2,K-1),(N-3,K-1)...(N-X,K-1) , (N-1,K-1) = (N-2,K-2),(N-3,K-2),(N-4,K-2)...(N-X,K-2) ,... (1,1) = (1) */
  CFactorials<double> * _factn ; /*! This will be created in constructor and used to get calculate the binomial coefficient of N,K pairs (the number of possible selections of K in size from set of N. */
  
	size_t _N ;  /*! The number in the set to choose from */
  size_t _K ;  /*! the number of values to choose from the the _N possiblities */
  size_t _max_steps ;  /*! = N - K + 1   - the extra one is to hold 0 as the first value in the scan. This is used to create vectors of the same size for each map input. Useful for a vectorised version of a scan of values. */
  

	/*! swaps N and K if K is greater than N 
	 *	 (not currently used)
	 */
  size_t MakeN_LTE_K (size_t N_IN)
  {
     size_t retVal ;
     if (N_IN < (this->_K-1))
        retVal = N_IN ;
      else 
        retVal = this->_K-1 ;
        
    return retVal ;
  }
  
  /*! 
   * \details Uses the recursive nature of the binomial theorm to compute the 'scan' of binomial sequence
   *  N.B. The scan of a set of numbers is the incremental sum of all the previous numbers in the set.
   * Example (N,K) = (8,4)  = sum [0, (7,3),(6,3),(5,3),(4,3),(3,3)]
   *                  70     = sum [0, 35,    20,   10,    4,   1   ]
   *                  result =     [0, 35,    55,   65,    69,  70  ]  = scan
   */
  std::vector<size_t> * ComputeNKVectorScan (size_t n_IN, size_t k_IN)
  {
    std::vector<size_t> * v_temp_scan = new std::vector<size_t>() ;
    v_temp_scan->reserve(_max_steps) ;
    
    v_temp_scan->push_back(0) ;  // first value is 0 in the scan results vector
    
    for (size_t t1 = n_IN ; t1 >= k_IN; t1--)
    {
     // if (k_IN == 1) std::cout << " (n,k) = " << t1<<  "," << k_IN << std::endl ; 
      size_t res = this->_factn->nchoosek(t1-1,k_IN-1)  ;  // -1 as it is recursive 
      v_temp_scan->push_back(res + v_temp_scan->back()) ;
    }
    
    // Fill the remaining scan vector values with the same value as the last computed up to max_steps
    size_t t_start_repeat = n_IN - k_IN + 1 ;
    for (size_t t1 = t_start_repeat; t1 < _max_steps-1; t1++)
    { 
      v_temp_scan->push_back(v_temp_scan->at(t1)) ;
    }
    
    return v_temp_scan ;
  }
  
  
  size_t ClosestDivisibleByNumThreads(size_t arraySizeIN, size_t arrayAlignemet)
  {
    double result = 0 ;
    result = (double) arraySizeIN / (double) arrayAlignemet ;
    result = ceil(result) ;
    return (size_t) ( result * arrayAlignemet ) ;
  }
    
  
public: 
  /*!
   * Creator - 
	 * \param nk_pairIN the binomial (N,K) pair of numbers - N being the Number of choices and K being how many to select from the N.
	 * \param arrayAlignemetIN - used to specify what multiple an arry should be extended to
   */
  CMapSelectKFromN( std::pair<size_t, size_t> nk_pairIN, size_t arrayAlignemetIN ) 
  {
    
    this->_factn = (CFactorials<double> *) new CFactorials<double>(20) ;
    
    this->_N = nk_pairIN.first ;
    this->_K = nk_pairIN.second ;
    
	//if (this->_N < this->_K)
   //  _PRINTSTD << "CMapSelectKFromN() Creator N=" << this->_N << " K=" << this->_K << ") " << std::endl  ;
  
    assert(this->_N >= this->_K) ;
    assert(this->_N >= 0) ;
    assert(this->_K >= 0) ;

    // The longest vector will be (N-1,k-1) - e.g. if N = 8, K=4 : (7,3) = 0, 35, 55, 65, 69, 70, so _arraySize = 6  including the +1 for the initial 0 at the start
    size_t _arraySize = this->_N - this->_K + 2;
    this->_max_steps = ClosestDivisibleByNumThreads(_arraySize, arrayAlignemetIN ) ;
    
    
    key_NK_pair t_pair ;  // this is a std::pair<size_t, size_t>
    
    size_t sub_step = 0 ;
    for (size_t K_step = this->_K; K_step >= 1 ; K_step-- )    
    {
      for (size_t N_step = this->_N-sub_step ; N_step  >= K_step ; N_step-- )
      { 
        t_pair = std::make_pair(N_step,K_step) ;
        std::vector<size_t> * v_nk = this->ComputeNKVectorScan(N_step, K_step) ;
        this->map_nchoosek.insert(std::make_pair(t_pair, v_nk)); // creates the map to pointers of the scan vectors
      }
      sub_step++ ;
    } 

	 t_pair = std::make_pair(0,0) ;
	 std::vector<size_t> * v_nk = new std::vector<size_t>(this->_max_steps,1) ;
	 size_t * ptr_data = v_nk->data()  ;
	 ptr_data[0] = 0 ;
    // std::vector<size_t> * v_nk = this->ComputeNKVectorScan(0, 0) ;
     this->map_nchoosek.insert(std::make_pair(t_pair, v_nk)); // creates the map to pointers of the scan vectors


  } // end of creator
  
  /*!
   * destructor
   */ 
   ~CMapSelectKFromN( void ) 
   {    
      key_NK_pair t_pair ;  // this is a std::pair<size_t, size_t>
      
      // iterate over the map and release the memory for each vector<> scan data present
      size_t sub_step = 0 ;
      for (size_t K_step = this->_K; K_step >= 1 ; K_step-- )
      {
        //sub_step++ ;
        for (size_t N_step = this->_N-sub_step ; N_step >= K_step ; N_step-- )
        {
          t_pair = std::make_pair(N_step,K_step) ;
          std::vector<size_t> * v_nk = this->map_nchoosek[t_pair] ;
          delete v_nk ;
        }
        sub_step++ ;
      }
      
  } // end of destructor
  
  /*! accessor function
   \return Returns the N member value, which represents the total number of possible choices in a set */
  size_t GetN() { return this->_N ; }
  /*! accessor function 
   * \return Returns the K member value, which represents the number selections wanted from a set */
  size_t GetK() { return this->_K ; }
  
  
  /*! Adds extra <N,K> pairs to the original map 
   * 
	 * \param nk_pairIN the binomial (N,K) pair of numbers - N being the Number of choices and K being how many to select from the N.
	 * \param arrayAlignemetIN - used to specify what multiple an arry should be extended to
   */
	void AddPairs( size_t nIN, size_t kIN, size_t arrayAlignemetIN ) 
  {
		std::pair<size_t, size_t> nk_pair ;
		nk_pair.first = nIN ;
		nk_pair.second = kIN ;
		this->AddPairs(nk_pair, arrayAlignemetIN) ;
	}

	 /*! This will add extra <N,K> pairs to the original map 
	 * \param nk_pairIN the binomial (N,K) pair of numbers - N being the Number of choices and K being how many to select from the N.
	 * \param arrayAlignemetIN - used to specify what multiple an arry should be extended to
   */
  void AddPairs( std::pair<size_t, size_t> nk_pairIN, size_t arrayAlignemetIN ) 
  { 
   // this->_factn = (CFactorials<double> *) new CFactorials<double>(20) ;
    this->_N = nk_pairIN.first ;
    this->_K = nk_pairIN.second ;
    
    assert(this->_N >= this->_K) ;
    assert(this->_N >= 0) ;
    assert(this->_K >= 0) ;
    
    size_t _arraySize = this->_N - this->_K + 2; // The longest vector will be (N-1,k-1) - e.g. if N = 8, K=4 : (7,3) = 0, 35, 55, 65, 69, 70  with + 1 for the initial 0 at the start.
    this->_max_steps = ClosestDivisibleByNumThreads(_arraySize, arrayAlignemetIN ) ;
    
    key_NK_pair t_pair ;  // this is a std::pair<size_t, size_t>
    
    size_t sub_step = 0 ;
    for (size_t K_step = this->_K; K_step >= 1 ; K_step-- )    
    {
      for (size_t N_step = this->_N-sub_step ; N_step  >= K_step ; N_step-- )
      { 
        t_pair = std::make_pair(N_step,K_step) ;
        std::vector<size_t> * v_nk = this->ComputeNKVectorScan(N_step, K_step) ;
        this->map_nchoosek.insert(std::make_pair(t_pair, v_nk)); // creates the map to pointers of the scan vectors
      }
      sub_step++ ;
    } 
  }
	

  /*!
   * \details This is an attempt to create a vectorisable version of the array index proceedure
   * Each of the input seq_numIN has a unique answer
   * \param seq_numIN is a 1 based request - generally will be a randomly generated number between 1 and the number of possible 
	 *	combinations as given by the binomial coeficient (N,K) (from N values, select K). The upper limit can be very large, so may 
	 *  require more bits than given by size_t
	 *				
	 * \return Return vector are the (K) array indexes to be choosen from an array of (N) values
   */ 
  std::vector<size_t> GetScan_Vectorised(size_t seq_numIN)
  {
	  assert(seq_numIN > 0) ;

    std::vector<size_t> res(this->_K, 0) ; // fill the result vector with 0's - we are selecting K at a time, so the result vector will be of size K
    std::vector<char>   resbool(this->_max_steps,0) ; // std::vector<bool> is actually a specialised STL container that uses a bitwise storage of [t | f] - may not be comutationally efficent.
    std::vector<size_t> * scan ; 
    key_NK_pair t_pair ; 
    
    
    size_t N_step = this->_N ;
    size_t K_step = this->_K ;
    size_t input = seq_numIN - 1 ;
    char * ptr_resbool  = NULL ;
    size_t * ptr_res = res.data() ;
    size_t add_extra = 0 ;
    size_t * ptr_scan = NULL ; 
    size_t first_true = 0; 
    ptr_resbool = resbool.data() ;
    
    
    for (size_t level = 0 ; level < this->_K ; level++ )
    {
     // std::cout << "START OF LOOP "<< level << ":  \tN_step: " << N_step << ",\tK_step: " << K_step << "\tinput: " << input ;
      
      t_pair      = std::make_pair(N_step,K_step) ;
      scan        = this->map_nchoosek[t_pair] ;
  	  if (scan == NULL) 
  	  {
  		  _PRINTERROR << "GetScan_Vectorised Error(): scan data was returned as NULL - the std::pair<N,K> (N=" << N_step << " K=" << K_step << ") does not exist in the std::map" << std::endl  ;
  		  break ;
  	  }
      ptr_scan    = scan->data() ;
      
     // this->PrintScan(N_step,K_step) ;
      
      for (size_t test = 0 ; test < this->_max_steps ; test++ )  // set the resbool vector value to true or false dependent on if the input is less that the value in the scan vector
      {
        if (input < ptr_scan[test])
          ptr_resbool[test] = 0 ;
        else
          ptr_resbool[test] = 1 ; // we could set first_true here rather than calculating it in the next loop
      }
         
      first_true = 0; // -1 as the first value is always 'false' as 0 < 0 is false

      for (size_t test = 0 ; test < this->_max_steps ; test++ )  // sum the values
      {
          first_true += ptr_resbool[test] ;
          ptr_resbool[test] = 0 ; // set as 0 as we go
      }
            
      first_true = first_true - 1 ; // -1 as the first value is always 'false' as 0 < 0 is false
      ptr_res[level] = level + first_true + add_extra ;
  	  N_step = N_step - first_true -1  ;
  	  K_step = K_step - 1 ; 
      input = input - ptr_scan[first_true] ; 
      add_extra += first_true  ;
    }
        
    return (res) ;
  }

  /*!
   * \details This is an attempt to create speed thing up
   * Each of the input seq_numIN has a unique answer
   * \param seq_numIN is a 1 based request - generally will be a randomly generated number between 1 and the number of possible 
	 *	combinations as given by the binomial coeficient (N,K) (from N values, select K). The upper limit can be very large, so may 
	 *  require more bits than given by size_t
	 *				
	 * \return Return vector are the (K) array indexes to be choosen from an array of (N) values
   */ 
  std::vector<size_t> GetScan_shortcut(size_t seq_numIN)
  {
		assert(seq_numIN > 0) ;

    std::vector<size_t> res(this->_K, 0) ; // fill the result vector with 0's - we are selecting K at a time, so the result vector will be of size K
    std::vector<size_t> * scan ; 
    key_NK_pair t_pair ; 
       
    size_t N_step = this->_N ;
    size_t K_step = this->_K ;
    size_t input = seq_numIN - 1 ;
    size_t * ptr_res = res.data() ;
    size_t add_extra = 0 ;
    size_t * ptr_scan = NULL ; 
    size_t first_true = 0; 
   
    
    for (size_t level = 0 ; level < this->_K ; level++ )
    {

    //  _PRINTSTD << "START OF LOOP "<< level << ":  \tN_step: " << N_step << ",\tK_step: " << K_step << "\tinput: " << input << std::endl ;

      
      t_pair      = std::make_pair(N_step,K_step) ;
      scan        = this->map_nchoosek[t_pair] ;
      
  	  if (scan == NULL) 
  	  {
  		  //_PRINTERROR << "GetScan_shortcut Error(): scan data was returned as NULL - the std::pair<N,K> (N=" << N_step << " K=" << K_step << ") does not exist in the std::map" << std::endl ;
  		  break ;
  	  }
      ptr_scan    = scan->data() ;   
     
      first_true = 1;
      size_t test ;
      test = 1 ;
      while ( input >= ptr_scan[test] )
      {
        first_true++ ;
		    test++ ;
      }
            
      first_true = first_true - 1 ; // -1 as the first value is always 'false' as 0 < 0 is false
      ptr_res[level] = level + first_true + add_extra ;

      //N_step = N_step - first_true -1  ;
      //K_step = K_step - 1 ;
	   N_step = N_step - first_true -1  ;
      if(  ((int)N_step - (int)first_true - 1 ) >= 0){
        N_step = N_step - first_true -1  ;
      } else {
         // Rcpp::warning("N_step is less than zero");
      }
        
      if(  (((int)K_step) - 1) >= 0){
        K_step = K_step - 1 ;
      } else {
       // Rcpp::warning("K_step is less than zero");
      }
      
	  
      input = input - ptr_scan[first_true] ; 
      add_extra += first_true  ;
    }
        
    return (res) ;
  }  // end of GetScan_shortcut() 

  /*!
   * \details This is an attempt to create speed thing up
   * Each of the input seq_numIN has a unique answer
   * \param seq_numIN is a 1 based request - generally will be a randomly generated number between 1 and the number of possible 
   *	combinations as given by the binomial coeficient (N,K) (from N values, select K). The upper limit can be very large, so may 
   *  require more bits than given by size_t
   *				
   * \return Return vector are the (K) array indexes to be choosen from an array of (N) values
   */ 
  std::vector<size_t> GetScan_shortcut_debug(size_t seq_numIN, std::string str_fromwhere)
  {
    assert(seq_numIN > 0) ;
    
    std::vector<size_t> res(this->_K, 0) ; // fill the result vector with 0's - we are selecting K at a time, so the result vector will be of size K
    std::vector<size_t> * scan ; 
    key_NK_pair t_pair ; 
    
    size_t N_step = this->_N ;
    size_t K_step = this->_K ;
    size_t input = seq_numIN - 1 ;
    size_t * ptr_res = res.data() ;
    size_t add_extra = 0 ;
    size_t * ptr_scan = NULL ; 
    size_t first_true = 0; 
    size_t st_last ;
    
    for (size_t level = 0 ; level < this->_K ; level++ )
    {
      
      //  _PRINTSTD << "START OF LOOP "<< level << ":  \tN_step: " << N_step << ",\tK_step: " << K_step << "\tinput: " << input << std::endl ;
      
      assert(N_step >= K_step) ;
      t_pair      = std::make_pair(N_step,K_step) ;
      scan        = this->map_nchoosek[t_pair] ;
      
      if (scan == NULL) 
      {
        _PRINTERROR << "GetScan_shortcut_debug Error(): scan data was returned as NULL - the std::pair<N,K> (N=" << N_step << " K=" << K_step << ") does not exist in the std::map" << std::endl ;
        break ;
      }
      ptr_scan    = scan->data() ;   
      
      first_true = 1;

			st_last = 0 ;
       while (( input >= ptr_scan[first_true] ) && (first_true < (N_step-K_step+1)) )
      {
        first_true++ ;
      }
      
      first_true = first_true - 1 ; // -1 as the first value is always 'false' as 0 < 0 is false
      ptr_res[level] = level + first_true + add_extra ;
      
      if(  ((long long int)N_step - (long long int)first_true - 1 ) >= 0){
        N_step = N_step - first_true -1  ;
      } else {
        _PRINTERROR << "GetScan_shortcut_debug() Error: seqnum: " << seq_numIN << ", " << str_fromwhere << ". N_step is less than zero: level=" << level << " N=" << N_step  << " first_true=" << first_true  <<  " K=" << K_step  << std::endl ; 
       // Rcpp::warning("N_step is less than zero" );
      }
      
      if(  (((long long int)K_step) - 1) >= 0){
        K_step = K_step - 1 ;
      } else {
        _PRINTERROR << "GetScan_shortcut_debug() Error: " << seq_numIN << ", " << str_fromwhere << ". K_step is less than zero: level=" << level << " N=" << N_step  << std::endl ;
      }
      
      
      input = input - ptr_scan[first_true] ; 
      add_extra += first_true  ;
    }
    
    return (res) ;
  }  // end of GetScan_shortcut_debug() 
  
  
  
  
  double nchoosek(size_t nIN, size_t kIN) 
  {
    return (this->_factn->nchoosek(nIN,  kIN)) ;
  }
  
  
  
  void PrintScan(size_t N_stepIN, size_t K_stepIN )
  {
    key_NK_pair t_pair ;
    t_pair = std::make_pair(N_stepIN,K_stepIN) ;
    std::vector<size_t> * v_nk = this->map_nchoosek[t_pair] ;
    std::cout << " ( " << N_stepIN << ",\t" << K_stepIN << " ) :\t "  ;
    for (size_t t1 = 0 ; t1 < v_nk->size() ; t1++)     
    {
      std::cout << v_nk->at(t1) << "\t" ;
    }
    std::cout << std::endl ;
  }
  
  // Print vectors of scan data present in the std::map container
  void PrintScans(size_t n_IN, size_t k_IN)
  {
    PrintScans( n_IN, k_IN, true) ;
  }
  
  /*! Print vectors of scan data present in the std::map container that belong to the (N,K) binomial series
	* N.B. there may be more than one set of binomial series (N,K) sets in the map and the sets may have
	* shared items (i.e. the map contains the union of a number of (N,K) binomial series)
	*/
  void PrintScans(size_t n_IN, size_t k_IN, bool printAll)
  {
    assert(n_IN >= k_IN) ;
    assert(k_IN <= this->_K) ;
    assert(n_IN <= this->_N) ;
    
    if (printAll == true)
    {
      size_t sub_step = 0 ;
      for (size_t K_step = this->_K; K_step >= 1 ; K_step-- )
      {
        for (size_t N_step = this->_N-sub_step ; N_step >= K_step ; N_step-- )
        {
          PrintScan(N_step,K_step);
        }
        sub_step++ ;
      } 
    }
    else // just print the scan of interest
    {
      PrintScan(n_IN,k_IN) ;
    }
  }

	/*! Prints all vectors of scan data present in the std::map container
	*  The order of priniting is not particulalry clear to be at the moment, but conforms to what the ::iterator returns
	*/
  void PrintAllScans( void )
  {   
		 key_NK_pair t_pair ;
		 for (std::map<key_NK_pair, std::vector<size_t >*>::iterator it=this->map_nchoosek.begin(); it!=this->map_nchoosek.end(); ++it)
		 {
			  t_pair = it->first ; 
				std::vector<size_t> * v_nk = it->second ;
				std::cout << " ( " << t_pair.first << ",\t" << t_pair.second << " ) :\t "  ;
				for (size_t t1 = 0 ; t1 < v_nk->size() ; t1++)     
				{
					std::cout << v_nk->at(t1) << "\t" ;
				}
				std::cout << std::endl ;
		 }

  }

} ; // end of class CMapSelectKFromN
