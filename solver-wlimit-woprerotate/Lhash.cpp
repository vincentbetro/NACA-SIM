#include "Mesh.h"

void Mesh::Lhash()
  {
  int i, j, k, flag; //counters

  //allocate for hash the first time through
  if (hash == 0)
    {
    hash = new List*[nn]; //each node has a list of nodes it's attached to
    for (i = 0; i < nn; i++)
      hash[i] = new List(); //make each node's list
    //printf("\nAllocating for hash.");
    }
  else
    {
    //we just redimension each row to zero (since nn still the same, just different dependencies now)
    for (i = 0; i < nn; i++)
      hash[i]->Redimension(0);
    //printf("\nOnly redimensioning hash.");
    }

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

  //now, determine size of mat
  mdim = 0; //number of entries in hash table
  for (i = 0; i < nn; i++)
    {
    mdim+=hash[i]->max;
    }

  return;
  }
