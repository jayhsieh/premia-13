#ifndef _BSND_PATH_H
#define _BSND_PATH_H

typedef struct PremiaBSnd
{
  PnlVect *spot; /*!< spot vector */
  PnlMat *LGamma; /*!< Cholesky factorization of the correl */
  PnlVect *sigma; /*!< vector of volatilities */
  double r; /*!< interest rate */
  PnlVect *divid; /*!< dividend rate */
  int d; /*!< size of the model */
  int timesteps; /*!< number of timesteps to use */
} PremiaBSnd;

extern void premia_bs_path (PnlMat *path, const PremiaBSnd *mod ,
                            const PnlMat *G, double T, const PnlVect *drift);

#endif /* _BSND_PATH_H */
