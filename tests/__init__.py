# -*- coding: utf-8 -*-
"""Unit test package for gffpandas."""

# standard library imports
import functools


def print_docstring():
    """Decorator to print a docstring."""

    def decorator(func):
        """Define decorator"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Print docstring and call function"""
            print(func.__doc__)
            return func(*args, **kwargs)

        return wrapper

    return decorator
