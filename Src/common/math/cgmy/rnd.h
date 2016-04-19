#ifndef _RND_H
#define _RND_H

class StableRnd{
  float B;
  float S;
  float alpha; 
  float sigma; 
  float beta; 
  float mu;
  int generator;
 public:
  StableRnd(float alpha, float sigma, float beta, float mu,int generator);
  float next();
};

float stablernd(float alpha, float sigma, float beta, float mu,int generator);

#endif
