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

class LMM
/* this class is the only one to contain information about the LIBOR's rates. We may value other bermudean derivative by changing only this class */
{
 private:  
  double sigma;//constant volatility
  double delta;
  double theta;//strike
  int p;//precision of the log-Euler scheme
  double dt, sqrt_dt; //time space (and its squared root) for the log-Euler scheme
  double L_0;//flat initial value for the rates 
  valarray<double> actual; //actualization
  valarray<double> B; //zero coupon bonds

 public :
  double payoff;
  int first; //beginning of the swaption
  int k; //number of remaining exercise dates   
  int kinit; //initial value of k
  valarray<double> L; //LIBOR rates

 public:
  LMM(int _first, int _k, double _L_0, double _sigma, double _delta, int _p, double _theta);

  void Scheme();
  void InitialCond();
  
  double CF();
  double next_euro();
  void print_MC(int N);
  void print_CF();

  double trajup();
  void estimateup(int N);
};

class Y1 : public LMM
{ 
 private :
  int N;

 public:
  Y1(int first, int k, double L_0,  double sigma, double delta, int p, double theta, int _N);

  double short_traj1();
  void estimate();
};

class Y2 : public LMM
{  
 private :
  int N1;
  int N2;
  valarray<double> P; 
  valarray<double> Lmem;

 public:
  Y2(int first, int _k, double L_0, double sigma, double delta, int p, double theta, int _N1, int _N2);

  void traj1();
  double best_y1();
  double short_trajZ();
  void estimate();//Y2 - Y1
};

class Andersen : public LMM
{
 private :
  int strat;
  int div;
  int NH;
  int NS;
  valarray<valarray<double> > po;
  valarray<valarray<double> > euro;// for the strategies 2, 3, 4 and 5

 public :
  valarray<double> H;

 public :
  Andersen(int _strat, int first, int _k, double L_0,  double sigma, double delta, int p, double theta, int _div, int _NH, int _NS);

  void tab_payoff(int N);
  void Move_H(int step, double h);
  double golden(double ax, double tol, int step);
  double Evaluate_H();  
  void estimate();
  void Find_H();
};

class Andersenhat : public Andersen
{
 private :
  int N1;
  int N2;
  valarray<double> P; 
  valarray<double> Lmem;//memory

 public :
  Andersenhat(int _N1, int _N2, int strat, int first, int k, double L_0, double sigma, double delta, int p, double theta, int div, int NH);

  void trajA();
  double best_yA();
  double short_trajZ();
  void estimate2();
};

class Y1up: public LMM
{
 private :
  int NK;
  int N1up;
  int N2up;
  valarray<double> Lmem1;//memory
  valarray<double> Lmem2;//memory
  valarray<double> P; 

  double up;
  double low;

 public :
  Y1up(int first, int k, double L_0, double sigma, double delta, int p, double theta, int _NK, int _N1up, int _N2up);
  
  double short_traj1();
  double e_v_short_traj1();
  double expected_value();
  double traj1up();
  double alpha_determination();
  void estimate();
};

class YAup : public Andersen
{
 private :
  int NK;
  int N1up;
  int N2up;
  valarray<double> Lmem1;//memory
  valarray<double> Lmem2;//memory
  valarray<double> P; 

  double up;
  double low;

 public :
  YAup(int strat, int first, int k, double L_0, double sigma, double delta, int p, double theta, int div, int NH, int _NK, int _N1up, int _N2up);

  double short_trajA();
  double e_v_short_trajA();
  double expected_value();
  void trajAup();
  double alpha_determination();
  void estimate2();
};
