#!/usr/bin/awk -f

BEGIN{ }
{  
            if( NR%2 == 0 ) print $0 
}
