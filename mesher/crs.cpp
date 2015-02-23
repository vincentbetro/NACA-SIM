#include "crs.h"

// hash takes number of nodes (dimension of List hash), number of elements in hash, hash table, and blank ia, ja, iau to be passed out full

void crs(int nn, int mdim, int *ia, int *ja, int *iau, List **hash)
  {
  int i, j, k; //counters

  //first, we init ia, based on the dimensions of each row of our hash table
  ia[0] = 0; //node 0 's connectivity always begins with the 0th element in ja
  for (i = 1; i < nn; i++)
    {
    ia[i] += ia[i-1] + hash[i-1]->max; //we add the previous index with the number of elements in the previous row to get the starting index for the next row
    }
  ia[nn] = mdim; //the last element in ia is the number of elements in ja

  //now, init ja by using the hash table, loop thru each row in hash table
  for (i = 0; i < nn; i++)
    {
    k = 0; //start hash table column counter over at 0 each time
    //loop through all entries of ja, indexing off ia, and insert node i 's entries
    for (j = ia[i]; j < ia[i+1]; j++)
      {
      ja[j] = hash[i]->list[k];
      k++; //increment k, and since ia is based off the entries in hash, it will only let k increment as many entries as there are in hash
      }
    }

  //finally, init iau by looking for each node in order, taking them from their "section" of ja, thus these are the diagonal elements in mat
  for (i = 0; i < nn; i++)
    {
    for (j = ia[i]; j < ia[i+1]; j++)
      {
      if (ja[j] == i)
        {
        iau[i] = j;
        }
      }
    }

  /*//now, print these for debugging
  printf("\nia:\n");
  for (i = 0; i < nn+1; i++)
    {
    printf("%d ",ia[i]);
    }
  printf("\n");

  printf("\nja:\n");
  for (i = 0; i < mdim; i++)
    {
    printf("%d ",ja[i]);
    }
  printf("\n");

  printf("\niau:\n");
  for (i = 0; i < nn; i++)
    {
    printf("%d ",iau[i]);
    }
  printf("\n");*/

  /*//debug
  List checker;
  for (i = 0; i < nn; i++)
    {
    checker.Check_List(iau[i]);
    }
  printf("\nchecker.max = %d, nn = %d\n",checker.max,nn);
  for (i = 0; i < nn; i++)
    {
    if(ja[iau[i]] != i)
      printf("\nYou have an index discrepancy at node %d\n",i);
    }
  for (i =0; i < nn; i++)
    {
    for (j = ia[i]; j < ia[i+1]; j++)
      {
      if (!hash[i]->Is_In_List(ja[j]))
        {
        printf("\nYou are missing node %d from list %d!\n",ja[j],i);
        }
      }
    }*/
        
  return;
  }
