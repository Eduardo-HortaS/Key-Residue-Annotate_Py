"""
decorators.py

Copyright 2024 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.

This script contains all decorator functions used in the SSF-Predict pipeline.

"""

import time
import logging
from typing import Callable, Any
import memory_profiler

def measure_time(func: Callable[..., Any], logger: logging.Logger):
    """
    Measures the time it took to execute the decorated function.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"\n{func.__name__} execution time: {execution_time:.6f} seconds\n")
        return result
    return wrapper

def measure_memory(func: Callable[..., Any], logger: logging.Logger):
    """
    Measures the memory usage before and after the decorated function's execution.
    """
    def wrapper(*args, **kwargs):
        memory_before = memory_profiler.memory_usage()[0]
        result = func(*args, **kwargs)
        memory_after = memory_profiler.memory_usage()[0]
        memory_used = memory_after - memory_before
        logger.info(f"\nMemory increase for {func.__name__}: {memory_used:.2f} MB\n")
        return result
    return wrapper

def measure_time_and_memory(func: Callable[..., Any], logger: logging.Logger):
    """
    Measures both execution time and memory usage increase of the decorated function.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        memory_before = memory_profiler.memory_usage()[0]
        result = func(*args, **kwargs)
        end_time = time.time()
        memory_after = memory_profiler.memory_usage()[0]
        execution_time = end_time - start_time
        memory_used = memory_after - memory_before
        logger.info(f"\n{func.__name__} execution time: {execution_time:.6f} seconds")
        logger.info(f"Memory increase for {func.__name__}: {memory_used:.2f} MB\n")
        return result
    return wrapper
