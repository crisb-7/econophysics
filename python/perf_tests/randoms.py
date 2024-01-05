import time
import random
import numpy as np

N = int(1e7)

def bench_time_loop(N):
    def decorator(func):
        def wrapper(*args, **kwargs):
            times = [None]*N
            for i in range(N):
                start_time = time.time()
                func(*args, **kwargs)
                elapsed_time = time.time() - start_time
                times[i] = elapsed_time
            print(f"{func.__name__} time (mean):", np.mean(times), "seconds")
        
        return wrapper
    return decorator


# Simple random, Numpy vs native Python
@bench_time_loop(N)
def numpy_random():
    numpy_rand = 2*np.random.rand()-1

@bench_time_loop(N)
def python_random():
    py_rand = 2*random.random() -1

def simple_random():
    python_random()
    numpy_random()
    # 10,000,000 iterations,
    # python_random time (mean): 2.564e-07 seconds
    # numpy_random time (mean): 5.570e-07 seconds

@bench_time_loop(N)
def random_choices_sign():
    random.choice([-1, 1])

@bench_time_loop(N)
def numpy_sign():
    np.sign(2*np.random.rand()-1)

def choices_sign():
    random_choices_sign()
    numpy_sign()
    # 10,000,000 iterations,
    # random_choices_sign time (mean): 5.6486e-07 seconds
    # numpy_sign time (mean): 1.3562e-06 seconds

if __name__ == "__main__":
    choices_sign()