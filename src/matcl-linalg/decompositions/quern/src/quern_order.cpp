#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#if 0
#include "quern.h"

int QUERN_get_rbfs_column_ordering(int m,
                                   int n,
                                   const int* A_row_start,
                                   const int* A_column_index,
                                   int* column_order)
{
   if(m<=0 || n<=0 || !A_row_start || !A_column_index || !column_order)
      return QUERN_INPUT_ERROR;
   // get some memory to work in
   int* work=(int*)std::malloc((n+1+A_row_start[m]+n+n)*sizeof(int));
   if(!work) return QUERN_OUT_OF_MEMORY;
   // figure out number of entries in each column
   int* column_start=work;
   std::memset(column_start, 0, (n+1)*sizeof(int));
   for(int i=0; i<m; ++i){
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j)
         ++column_start[A_column_index[j]+1];
   }
   // cumulative sum to get column_start
   for(int i=0; i<n; ++i)
      column_start[i+1]+=column_start[i];
   assert(column_start[n]==A_row_start[m]);
   // list the columns now
   int* row_index=column_start+(n+1);
   int* column_pointer=row_index+column_start[n];
   std::memcpy(column_pointer, column_start, n*sizeof(int));
   for(int i=0; i<m; ++i){
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j){
         int c=A_column_index[j];
         row_index[column_pointer[c]++]=i;
      }
   }
   // set up marker for BFS
   char* column_marker=(char*)(column_pointer+n);
   std::memset(column_marker, 0, n);
   // and do as many BFS as we need to hit all connected components
   int p=n;
   for(int root=0; root<n; ++root) if(!column_marker[root]){
      column_order[--p]=root;
      column_marker[root]=1;
      for(int i=p; i>=p; --i){
         int j=column_order[i];
         // add unmarked neighbour columns of j to ordering
         for(int k=column_start[j]; k<column_start[j+1]; ++k){
            int r=row_index[k];
            for(int a=A_row_start[r]; a<A_row_start[r+1]; ++a){
               int nbr=A_column_index[a];
               if(!column_marker[nbr]){
                  column_order[--p]=nbr;
                  column_marker[nbr]=1;
               }
            }
         }
      }
   }
   assert(p==0);
   std::free(work);
   return QUERN_OK;
}

#endif