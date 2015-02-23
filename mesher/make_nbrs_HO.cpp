#include "make_nbrs_HO.h"

//this routine takes a list of triangles to make nabors from, number of nodes, number of triangles, expanded tri array (even though we only use 0, 1, 2), nabor array (passed in init to -1 and out full), and node_node hash table
void make_nbrs_HO(Linked_List *mknbr, int nn, int nt, int **tri, int **nbr, Linked_List **Nhash)
{
  //declare linked node to move through mknbr and keep track of triangles that have been fully "naborized"
  Linked_Node *hd;
  //counters, placeholders for node indices
  int c, m, n, s;
  int n0, n1;

  //set linked_node to head of mknbr and do loop while we have elements who have not been naborized
  while (hd = mknbr->head)
  {
    c = hd->data;  //holds data from linked list

    //loop thru all sides of c looking for nabors...if none, leave uninit
    for (s=0; s < 3; s++)
    {
      //switch over nodes on each side of c to compare with nodes on sides of other triangles
      switch (s)
      {
        case 0: n0 = tri[c][0]; n1 = tri[c][1]; break;
        case 1: n0 = tri[c][1]; n1 = tri[c][2]; break;
        case 2: n0 = tri[c][2]; n1 = tri[c][0]; break;
      }
      //init m flag to -1 ... if left this way, we are on bd and if not, we fill with nabor number
      m = -1;
      //linked nodes to read thru linked list
      Linked_Node *hd0, *hd1;
      // each element in n0's list is a triangle with that node in it
      hd0 = Nhash[n0]->head;
      //while there are elements that have n0 in them and m not set, keep looking
      while (hd0 && m < 0)
      {
        //cell from n0 Nhash
        n = hd0->data;
        //if the triangle we found not equal to the one we are looking for (c) 
        if (n != c)
        {
          //now, look to other node on same side's Nhash to see if it is also in element n
          hd1 = Nhash[n1]->head;
          //while there are elements in n1's Nhash and m unset
          while (hd1 && m < 0)
          {
            //if the element in n1's hash matches the element in n0's hash, we have en element that shares a side
            //set m (nabor) equal to that element
            if (hd1->data == n)
              m = n;
            //move to the next element in n1's hash, even if unneeded
            hd1 = hd1->next;
          }
        }
        //move to the next element in n0's hash, even if unneeded
        hd0 = hd0->next;
      }
      //set triangle c's nabor on side s (based on switch) to m, even if still -1, meaning we are on bd
      nbr[c][s] = m;
    }
    //now that nabors have been established, remove c from list of triangles to be cycled thru
    mknbr->Remove(c);
  }

}
