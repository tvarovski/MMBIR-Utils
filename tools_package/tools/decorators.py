import functools

def count_calls(func):
    # decorator that counts the number of times a function is called

    #create a wrapper function that counts the number of times the function is called
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.num_calls += 1
        #print out the number of times the function has been called every 1000 times
        if wrapper.num_calls % 1000 == 0:
            print(f"{wrapper.num_calls} calls to {func.__name__}")
        return func(*args, **kwargs)
    
    #initialize the number of calls to 0
    wrapper.num_calls = 0

    return wrapper