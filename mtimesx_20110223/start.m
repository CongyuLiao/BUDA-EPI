mex -DDEFINEUNIX CFLAGS="\$CFLAGS -march=native" -largeArrayDims -lmwblas -lmwlapack -lgomp mtimesx.c 