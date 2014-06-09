/**
 * @Description Quicksort template function for an array of generic sortable items.        
* Operators <, >, and == must be defined for every type T in such a way that  
* (1) P((a<b) || (a>b) || (a==b)) = 1                                       
* (2) P( (a<b) && (a>b) ) = 0                                                 
* (3) P( (a<b) && (a==b) ) = 0                                                
* (4) P( (a>b) && (a==b) ) = 0                                                
* (5) a<b is equivalent to b>a                                                
* (6) ( (a<b) & (b<c) ) -> (a<c)                                              
* for every (a,b) in T*T                                                      
* and every (a,b,c) in T*T*T                                                  
*                                                                             
* Operator = must be, for every type T, a memberwise copy operator.           
*/

#ifndef SORT_ALGOS_H
#define SORT_ALGOS_H


namespace SortAlgos {
    
    typedef unsigned int Idx;
    template<class T> Idx qsPartition(T* items, const Idx start, const Idx end);
    template<class T> void quicksort(T* items, const Idx start, const Idx end);
    template<class T> void swap(T* item1, T* item2);
}


/**
 * @description sorts an array of items using quicksort algorithm.
*                  After calling this function, if i<j then the i-th item is
*                  < or == to the j-th item.
*  @param  pointer to array of sortable items(T*), lower array index(Idx), upper array index(Idx).
 * @return changes the object internally (void)
*/
template<class T> void SortAlgos::quicksort(T* items, const Idx start, 
                                            const Idx end) {
    if(start < end) {
        Idx j = qsPartition(items, start, end);
        quicksort(items, start, j);
        quicksort(items, j+1, end);
    }
}

/**
 * @description  partitions the items between 'start' and 'end' (extremes
*                  included) into two subarrays such that all items in the
*                  first subarray, starting at index 'start', are < or =
*                  to all items in the second subarray, ending at index 'end'.
*@param pointer to array of sortable items(T*), lower array index(Idx), upper array index(Idx).
 * @return the new array index(Idx)
 */
template<class T> SortAlgos::Idx SortAlgos::qsPartition(T* items, 
                                            const Idx start, const Idx end) {
    const T pivot = items[start];
    Idx i = start-1;
    Idx j = end+1;
    for(;;) {
        do --j; 
        while(items[j] > pivot);
        do 
            ++i;
        while(items[i] < pivot);
        if(i<j)
            swap(&items[i], &items[j]);
        else
            return j;
    }    
}


/**
* @description  monolithic swap between two items of type T.
*@param  pointers to each array of items(T*).
*@return changes the object internally (void)
*/
template<class T> void SortAlgos::swap(T* item1, T* item2) {
    T temp = *item1;
    *item1 = *item2;    
    *item2 = temp;
}

#endif
