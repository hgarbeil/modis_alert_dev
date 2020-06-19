
#include "median.h"

// a is the input data array and n is the number of points in the array										  
                                                               
unsigned short tir_median (unsigned short *a, int n )                                                     
{                                                           
       int i,j;                                              
       unsigned short median,t;                                    
                                                               

for (i = 1 ; i <= n-1 ; i++)                            
    {     /* Trip-i begins  */
	for (j = 1 ; j <= n-i ; j++)                         
    {                                                    
        if (a[j] <= a[j+1])                              
        { /* Interchanging values */
                     
            t = a[j];                                      
            a[j] = a[j+1];                                 
            a[j+1] = t;                                    
        }                                                
        else
                    continue ;
    }                                                    
    } /* sorting ends */
	/* calculation of median  */
	if ( n % 2 == 0)                                        
          median = (a[n/2] + a[n/2+1])/2.0 ;                   
    else                                                    
          median = a[n/2 + 1];                                 
                                                               
    return median ;                                                    
   }    
   
   
float fmedian(float *a, int n )
{                                                           
       int i,j;                                              
       float median,t ;
       int medind;
       medind = (int) ((float)n * 0.5) ;
                                                               



    /* Sorting begins */
for (i = 1 ; i <= n-1 ; i++)                            
    {     /* Trip-i begins  */
	for (j = 1 ; j <= n-i ; j++)                         
    {                                                    
        if (a[j] <= a[j+1])                              
        { /* Interchanging values */
                     
            t = a[j];                                      
            a[j] = a[j+1];                                 
            a[j+1] = t;                                    
        }                                                
              else continue ;
    }                                                    
    } /* sorting ends */
	/* calculation of median  */
	if ( n % 2 == 0)                                        
          median = (a[medind] + a[medind+1])/2.0 ;
    else                                                    
          median = a[medind ];
                                                               
    return median ;


}    
    

float tir_fmedian_avg (float *a, int npts )
{
int i ;
float total=0. ;
for (i=0; i<npts; i++) {
    total += a[i] ;

}
return total/npts ;
          
}
