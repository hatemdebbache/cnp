"""
    Usfl Functions Module
    
    useful function that might be needed in future projects
    
    MADE BY: Hatem D.
"""

#==================================================
def counter(f, *args, **kwargs):
    """
    Evaluate the approx. time of running a function

    Parameters
    ----------
    f : function namespace
        Specify the function to be tested.
    *args : List of arguments to be passed to 'f' without keywords.
    **kwargs : List of keyworded arguments to be passed to 'f'

    Returns
    -------
    t : time delta between function call and return
    ret : the return of 'f'

    """
    from datetime import datetime
    
    t1 = datetime.now()
    ret = f(*args, **kwargs)
    t2 = datetime.now()
    
    return t2-t1, ret 
#==================================================
