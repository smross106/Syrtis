"""
Miscellaneous tool functions for use elsewhere

"""

from numbers import Number

def is_numeric(x, positive=False, not_negative=False, unit=False):
    """
    Check if input x is either numeric, or is a list of numeric
    """
    if isinstance(x, Number):
        if positive:
            # x must be greater than zero
            if x > 0:
                return(True)
            else:
                return(False)

        elif not_negative:
            # x must be greater than or equal to zero
            if x >= 0:
                return(True)
            else:
                return(False)

        elif unit:
            # x must be between 0 and 1
            if x>=0 and x<=1:
                return(True)
            else:
                return(False)

        else:
            return(True)

    elif type(x) == tuple or type(x) == list:
        all_items_numeric = all(is_numeric(item, positive, not_negative) for item in x)
        return(all_items_numeric)
    
    else:
        return(False)
