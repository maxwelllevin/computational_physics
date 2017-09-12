# -*- coding: utf-8 -*-
"""
Fibonacci Numbers
# Computational Physics

Created on Tue Sep 12 15:55:16 2017

@author: Maxwell
"""



def while_fibonacci(n, fib_list=[]):
    """ Returns a list containing 'n' fibonacci numbers. """
    # Start our fibonacci chain:
    fib_1, fib_2 = 1, 1
    
    # Declare our index counter
    index = 0
    
    # Loop while index (serves as a counter for #fib found) is < n
    while index < n:
        fib_list.append(fib_1) 
        fib_1 += fib_2
        index += 1
        
        # After incrementing index, check to ensure index < n
        if index < n:
            fib_list.append(fib_2)
            fib_2 += fib_1
            index += 1
    
    # Return our list of fibonacci numbers
    return fib_list



# ======= Test our code ======= #


fibs = while_fibonacci(10000)
print(fibs)