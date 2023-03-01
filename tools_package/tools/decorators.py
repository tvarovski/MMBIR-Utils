# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

import functools
import logging

logger = logging.getLogger('__main__.' + __name__)

def count_calls(func):
    # decorator that counts the number of times a function is called

    #create a wrapper function that counts the number of times the function is called
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.num_calls += 1
        #print out the number of times the function has been called every 1000 times
        if wrapper.num_calls % 1000 == 0:
            logging.info(f"{wrapper.num_calls} calls to {func.__name__}")
        return func(*args, **kwargs)
    
    #initialize the number of calls to 0
    wrapper.num_calls = 0

    return wrapper

def time_elapsed(func):
    import time
    # decorator that times the elapsed time of a function
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):

        #start the timer
        start = time.time()
        #call the function
        result = func(*args, **kwargs)
        #stop the timer
        end = time.time()
        #calculate the elapsed time
        elapsed = end - start
        #print the elapsed time rounded to 2 decimal places
        logging.info(f"\n{func.__name__} took {elapsed:.2f} seconds to run...")
        return result
    
    return wrapper

def fancy_status(func):
    # decorator that prints a fancy message when a function is called
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logging.info(f"\n{'*'*5}Calling {func.__name__}...{'*'*5}\n")
        result = func(*args, **kwargs)
        logging.info(f"\n{'*'*5}Done with {func.__name__}!{'*'*5}")
        return result
    
    return wrapper
