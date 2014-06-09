#!/usr/bin/awk -f

BEGIN{ }
{  
            if( NR%2 ) print $0 
}
