# -*- coding: utf-8 -*-
"""
Exercise 2.12: Prime Numbers
(The slow way.. (not a sieve))

Created on Sun Sep 17 18:47:25 2017

@author: Maxwell
"""

prime_list = [2]

def is_prime(test_num):
    """ Returns True if prime. False otherwise"""
    # If we've cached it, it is prime
    if test_num in prime_list: return True
    
    # If it is divisible by something we've cached and it is not in the cache, it is composite
    for prime in prime_list:
        if test_num % prime == 0: return False
    
    # This is unnecessary for testing done in order, necessary for random prime testing
    # If it is is divisible by a number up to the sqrt of itself, it is composite
    index = prime_list[-1]
    while index**2 < test_num + 1:
        if test_num % index == 0: return False
        index += 1
    
    # If it is none of these, it is prime 
    prime_list.append(test_num)
    return True
    


def find_primes_to_n (limit=1000):
    """ Returns a list of primes under the given limit"""
    # Test a numbers within our range
    for n in range(3, limit):
        is_prime(n)
        
    return prime_list


# ======= Test our function ======= #

import datetime as dt
before = dt.datetime.today()

primes = find_primes_to_n(10000)

after = dt.datetime.today()

print(primes)
print("Time elapsed:", after - before)