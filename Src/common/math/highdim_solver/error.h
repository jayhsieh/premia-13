
#define FDERROR(x) fprintf(stderr, "ERROR in " __FILE__ "(%d): " x, __LINE__);

#ifdef DEBUG
#  define FDDEBUG(args)                                       \
  do                                                          \
  {                                                           \
    printf("Debug message in " __FILE__ "(%d): ", __LINE__);  \
    printf args;                                              \
  }                                                           \
  while(0);
#else
#  define FDDEBUG(args)
#endif

#define ON_LASPACK_ERROR                                                      \
  if (LASResult() != LASOK)                                                   \
  {                                                                           \
    FDERROR("LASPack, ");                                                     \
    WriteLASErrDescr(stderr);                                                 \
  }                                                                           \
                                                                              \
  if (LASResult() != LASOK)

