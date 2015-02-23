#include <stdio.h>
#include <stdlib.h>

#ifndef List_h
#define List_h

class List
{
  public:
  List()
  { dim = max = 0; list = 0; }
  ~List()
  { dim = max = 0; if (list != 0) free(list); }
  // constructor called in general use and covered by voxel.construct()
  void print(FILE *outf)
  { 
    fprintf(outf,"\nInteger list dimension =%d",dim);
    fprintf(outf,"\nInteger list maximum index =%d",max);
    for (int i=0; i < max; i++)
      fprintf(outf,"\nlist(%d)= %d",i,list[i]);
  }
  int Redimension(int size);
  int Is_In_List(int n);
  int Check_List(int n);
  int Add_To_List(int n);
  int Delete_From_List(int n);
  int Replace(int m, int n);
  int Index(int n); //allows to find where in list..new 10/4/06
  
  int *list;
  int dim, max;
};

#endif
