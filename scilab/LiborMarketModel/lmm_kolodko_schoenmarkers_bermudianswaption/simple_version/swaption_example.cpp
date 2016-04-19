/*---------------------------------------*/
/*   BERMUDAN SWAPTION PRICER            */
/*    Kolodko & Schoenmakers Algorithm   */
/*    Andersen Algorithm                 */
/*    Tight upper estimator              */
/*                                       */
/*---------------------------------------*/
/*  Julien Bourgouint, Premia July 2006  */
/*  julien_bourgouint@yahoo.fr           */ 
/*---------------------------------------*/

#include "swaption.cpp"

int main()
{
  srand48 ((unsigned) time (NULL));
  
  double delta = .5; // tenor equal to 6 months
  double L_0 = .06;// initial flat value
  double theta = .06; // strike
  float swaptionMat= 1.0;        // (years)
  float swapMat= 4.0;            // (years)
  double sigma = .2;// volatily
  int p = 1; // precision for the log_Euler Scheme (precision increases with p since dt = delta/p)
  
  int k = (int)(swapMat/delta + .5);// last exercise date 
  int first = (int)(swaptionMat/delta + .5);// first exercise date

  int choice;  
  cout << "Bermudan swaption's pricing by Julien Bourgouint, July 2006." << endl << "  0. upper & lower bound using fast approach" << endl << endl << "  Tight Lower bound" << endl << "  1. Andersen Algorithm" << endl << "  2. Kolodko & Schoenmakers Algorithm" << endl  << endl << "  Tight Upper bound" << endl << "  3. with Andersen Algorithm" << endl << "  4. with Kolodko & Schoenmakers Algorithm" << endl;

  cout << "Choose  "; cin >> choice;
  cout << endl;  
  switch(choice)
    {
    case 0:
      {
	LMM Lib(first, k, L_0, sigma, delta, p, theta);
	Lib.InitialCond();
	Lib.print_CF();
	
	int N;
	cout << "Number of Monte Carlo Iterations for the lower bond"<<endl;
	cin >> N;
	LMM Lib2(first, k, L_0,  sigma, delta, p, theta);
	Lib2.InitialCond();
	Lib2.print_MC(N); 
	
  	cout << endl << "# # # # # # " << endl << "You can get a graph of (Monte Carlo & Close Form) european prices by running scilab and writting" << endl << "	getf(~/MC.sci)" << endl << "	plot2d(MC())" << endl<< "and" << endl << "	getf(~/CF.sci)" << endl << "	plot2d(CF())" << endl << "# # # # # # " << endl << endl;

	int Nup;
	cout << "Number of Monte Carlo Iterations the upper bond"<<endl;
	cin >> Nup;
	LMM Libup(first, k, L_0, sigma, delta, p, theta);
	Libup.estimateup(Nup);
	
	break;
      }
    case 1: 
      {
	int strat;
	cout << "Choose a strategy (1, 2, 3, 4 or 5) "; cin >> strat;
	int NH;// for the computation of the exercise boundary
	cout << "Number of Monte Carlo Iterations for the computation of the exercise boundary"<<endl;
	cin >> NH;

	int NS ;// for the computation of the corresponding price
	cout << "Number of Monte Carlo Iterations for the computation of the corresponding price"<<endl;
	cin >> NS;
	int div = 7; // H is piecewise affine and div is the number of pieces
	// In case div = k - first + 1, H is not constructed by interpolation and precision is optimal 
	Andersen Ander(strat, first, k, L_0, sigma, delta, p, theta, div, NH, NS);
	Ander.estimate();
	break;
      }  
    case 2:
      {
	
	
	int N;// for the price estimator after 1 iteration
	cout << "Number of Monte Carlo Iterations for the price estimator after 1 iteration"<<endl;
	cin >> N;
	Y1 y1(first, k, L_0, sigma, delta, p, theta, N);
	y1.estimate();

	int N1;// for the price estimator after 2 iterations
	cout << "Number of Monte Carlo Iterations for the price estimator after 2 iteration (less than in the first iteration)"<<endl;
	cin >> N1;

	int N2 ;// for inner MC estimations
	cout << "Number of Monte Carlo Iterations for inner MC estimations"<<endl;
	cin >> N2;
	Y2 y2(first, k, L_0, sigma, delta, p, theta, N1, N2);
	y2.estimate();
	break;
      }
    case 3 :
      {
	int strat = 1;
	int NH;// for the computation of the exercise boundary
	cout << "Number of Monte Carlo Iterations for the computation of the exercise boundary"<<endl;
	cin >> NH;
	int div = 7;// H is piecewise affine and div is the number of pieces. 

	cout << "Number of Monte Carlo Iterations is the product of the 3 following MC parameters"<<endl;

	int N1upA;// Number of Monte Carlo Iterations for the upper price estimator
	cout << "Number of Monte Carlo Iterations for the upper price estimator (ex.1000))"<<endl;
	cin >> N1upA;

	int N2upA ;// for inner MC estimations
	cout << "Number of Monte Carlo Iterations for inner MC (ex.100)"<<endl;
	cin >> N2upA;

	int NKA ;// for inner inner MC estimations
	cout << "Number of Monte Carlo Iterations for inner inner MC estimations (ex.10)"<<endl;
	cin >> NKA;

  
	YAup yaup(strat, first, k, L_0, sigma, delta, p, theta, div, NH, NKA, N1upA, N2upA);
	yaup.estimate2();
	break;
      }   
    case 4 :
      {
	 cout << "Number of Monte Carlo Iterations is the product of the 3 following MC parameters"<<endl;

	int N1up;// Number of Monte Carlo Iterations for the upper price estimator
	cout << "Number of Monte Carlo Iterations for the upper price estimator (ex.1000))"<<endl;
	cin >> N1up;

	int N2up ;// for inner MC estimations
	cout << "Number of Monte Carlo Iterations for inner MC (ex.100)"<<endl;
	cin >> N2up;

	int NK ;// for inner inner MC estimations
	cout << "Number of Monte Carlo Iterations for inner inner MC estimations (ex.10)"<<endl;
	cin >> NK;
	
	Y1up y1up(first, k, L_0, sigma, delta, p, theta, NK, N1up, N2up);
	y1up.estimate();
	break;
      }
    default: 
      {
	break;
      }  
    }
}
