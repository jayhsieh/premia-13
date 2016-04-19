#include "stack-c.h"
#include "../demo/price_cdo.h"

#define	MAX(a, b)	( (a < b) ? b : a )

int 			sci_cdo(char			*fname) 
{
    static int 		p_ncomp, m_ncomp, n_ncomp;
    static int 		p_nominal, m_nominal, n_nominal;
    static int 		p_dates, m_dates, n_dates;
    static int 		p_tranches, m_tranches, n_tranches;
    static int 		p_intensity, m_intensity, n_intensity;
    static int 		p_xrates, m_xrates, n_xrates;
    static int 		p_yrates, m_yrates, n_yrates;
    static int 		p_trecovery, m_trecovery, n_trecovery;
    static int 		p_recovery, m_recovery, n_recovery;
    static int 		p_tcopula, m_tcopula, n_tcopula;
    static int 		p_pcopula, m_pcopula, n_pcopula;
    static int		p_tmethod, m_tmethod, n_tmethod, tmethod;
    static int		p_pmethod, m_pmethod, n_pmethod;
    static int		p_price, p_dl, p_pl, m, n, size_price;
    static int 		minlhs=1, maxlhs=3;
    static int		minrhs=1, maxrhs=13;

    CheckRhs(minrhs,maxrhs);
    CheckLhs(minlhs,maxlhs);
    GetRhsVar(1, "i", &m_ncomp, &n_ncomp, &p_ncomp);
    GetRhsVar(2, "d", &m_nominal, &n_nominal, &p_nominal);
    GetRhsVar(3, "d", &m_dates, &n_dates, &p_dates);
    GetRhsVar(4, "d", &m_tranches, &n_tranches, &p_tranches);
    GetRhsVar(5, "d", &m_intensity, &n_intensity, &p_intensity);
    GetRhsVar(6, "d", &m_xrates, &n_xrates, &p_xrates);
    GetRhsVar(7, "d", &m_yrates, &n_yrates, &p_yrates);
    GetRhsVar(8, "i", &m_trecovery, &n_trecovery, &p_trecovery);
    GetRhsVar(9, "d", &m_recovery, &n_recovery, &p_recovery);
    GetRhsVar(10, "i", &m_tcopula, &n_tcopula, &p_tcopula);
    GetRhsVar(11, "d", &m_pcopula, &n_pcopula, &p_pcopula);
    GetRhsVar(12, "i", &m_tmethod, &n_tmethod, &p_tmethod);
    GetRhsVar(13, "i", &m_pmethod, &n_pmethod, &p_pmethod);
    m = MAX(m_tranches, n_tranches)-1;
    n = 1;
    //if ((p_tmethod == 6)||(p_tmethod == 7)) size_price = 2*m;
    size_price = m;
    tmethod = *(istk(p_tmethod));
    if(tmethod >=6)  size_price = 2*m;
    if(tmethod>=8)     size_price=m;
    if(tmethod>=9)     size_price=m+5;  
    if(tmethod>=10)    size_price=m+50;
  //if (p_tmethod == 7) size_price = 2*m;    
    //else size_price = m;
    CreateVar(14, "d", &size_price, &n, &p_price);
    CreateVar(15, "d", &size_price, &n, &p_dl);
    CreateVar(16, "d", &size_price, &n, &p_pl);
    price_cdo(istk(p_ncomp), 
	       stk(p_nominal), 
	       MAX(m_dates, n_dates), 
	       stk(p_dates), 
	       MAX(m_tranches, n_tranches), 
	       stk(p_tranches), 
	       stk(p_intensity), 
	       MAX(m_xrates, n_xrates), 
	       stk(p_xrates), 
	       stk(p_yrates), 
	       istk(p_trecovery), 
	       stk(p_recovery), 
	       istk(p_tcopula), 
	       stk(p_pcopula), 
	       istk(p_tmethod), 
	       istk(p_pmethod), 
	       stk(p_price), 
	       stk(p_dl), 
	       stk(p_pl)); 
    LhsVar(1) = 14;
    LhsVar(2) = 15;
    LhsVar(3) = 16;

    return 0;
}
