#include "Lhash.h"

//takes number of nodes, triangles, tri conn, and empty node to node hash to fill in

void Lhash(int nn, int nt, int **tri, List **hash)
  {
  int i, j, k, flag; //counters

  //now, go through all nodes and look for those nodes they are connected to using tri
  for (i = 0; i < nn; i++)
    {
    //look thru tri to find node, if it is in a triangle, add other nodes to List using Check_List; be sure to add node itself as well
    for (j = 0; j < nt; j++)
      {
      flag = 0; //reset flag
      for (k = 0; k < 3 && !flag; k++)
        {
        if (tri[j][k] == i)
          flag = 1; //we are in this triangle, so grab all nodes 
        }
      if (flag)
        {
        for (k = 0; k < 3; k++)
          hash[i]->Check_List(tri[j][k]); //check to be sure we have no dups (always will dup node i)
        }  
      } //end for each node
    } //end looking thru nodes, hash table created

  /*//debugging...prints hash table
  printf("\nHash table:\n");
  for (i = 0; i < nn; i++)
    {
    for (j = 0;  j < hash[i]->max; j++)
      {
      printf("%d ",hash[i]->list[j]);
      }
    printf("\n");
    }*/

  return;
  }
