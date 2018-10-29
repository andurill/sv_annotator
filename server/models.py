class SV(object):
    '''Define SV object that holds all key information'''
    def __init__(self, svtype, bkp1, bkp2, genes, site1, site2, description, connection):
       self.__svtype = svtype
       self.__bkp1 = bkp1
       self.__bkp2 = bkp2
       self.__genes = genes
       self.__site1 = site1
       self.__site2 = site2
       self.__description = description
       self.__connection = connection
    
    def __setBkp():
        try:
            (self.__chr1, self.__coord1) = self.__bkp1.split(":")
        except ValueError as e:
            return e
        return


class Error(Exception):
    '''Base class for other exceptions'''
    pass
class MissingCytoBandError(Error):
    '''Raised when no cytoband was identified for a breakpoint'''
    def __init__(self):
        Exception.__init__(self, "No cytobands identified for the breakpoint.")
class MultipleCytoBandError(Error):
    '''Raised when multiple cytoband were identified for a breakpoint'''
    def __init__(self):
        Exception.__init__(self, "Multiple cytobands identified for the breakpoint.")
class CanonicalTranscriptError(Error):
    '''Raised when canonical transcript not found'''
    def __init__(self):
        Exception.__init__(self, "Cannot find a canonical transcript for the breakpoint.")
class MultipleCanonicalTxError(Error):
    '''Raised when multiple canonical transcript versions exist'''
    def __init__(self):
        Exception.__init__(self,"Multiple versions of canonical transcript exists.")

'''  
       In[90]: class Error(Exception):
    ...:    """Base class for other exceptions"""
    ...: pass
    ...:
    ...: class ValueTooSmallError(Error):
    ...:    """Raised when the input value is too small"""
    ...: pass
    ...:
    ...: class ValueTooLargeError(Error):
    ...:    """Raised when the input value is too large"""
    ...: pass
    ...:


In[91]: class Error(Exception):
    ...:    """Base class for other exceptions"""
    ...:    error_val = 4
    ...: pass
    ...:
    ...: class ValueTooSmallError(Error):
    ...:    """Raised when the input value is too small"""
    ...: pass
    ...:
    ...: class ValueTooLargeError(Error):
    ...:    """Raised when the input value is too large"""
    ...: pass
    ...:

In[92]: a = Error()

In[93]: a
Out[93]: __main__.Error()

In[94]: a.error_val
Out[94]: 4

In[95]: a = ValueTooSmallError

In[96]: a
Out[96]: __main__.ValueTooSmallError

In[97]: a = ValueTooSmallError()

In[98]: a
Out[98]: __main__.ValueTooSmallError()

In[99]: a = ValueTooSmallError

In[100]: a.error_val
Out[100]: 4

In[101]: a = ValueTooSmallError()

In[102]: a.error_val
Out[102]: 4

In[103]: a = 5

In[104]: try:
     ...: if a == 6:
     ...: raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
  File "<ipython-input-104-f1adcaaa94b2>", line 4
    except ValueTooLargeError as e:
         ^
SyntaxError: invalid syntax


In [105]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:         break
     ...:     except ValueTooLargeError as e:
  File "<ipython-input-105-c16b056201d9>", line 5
    except ValueTooLargeError as e:
         ^
SyntaxError: invalid syntax


In [106]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:         break
     ...:     except ValueTooLargeError:
  File "<ipython-input-106-0c2aec896aa0>", line 5
    except ValueTooLargeError:
         ^
SyntaxError: invalid syntax


In [107]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:         break
     ...: except ValueTooLargeError as e:
     ...:     print e
  File "<ipython-input-107-4e65e2ef3031>", line 4
    break
SyntaxError: 'break' outside loop


In [108]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?

In [109]: a
Out[109]: 5

In [110]: a=6

In [111]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [112]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooSmallError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?
---------------------------------------------------------------------------
ValueTooLargeError                        Traceback (most recent call last)
<ipython-input-112-2d1994a9d726> in <module>()
      1 try:
      2     if a==6:
----> 3         raise ValueTooLargeError
      4 except ValueTooSmallError as e:
      5     print e

ValueTooLargeError:

In [113]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [114]: a=5

In [115]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooSmallError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?

In [116]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?

In [117]: a=5

In [118]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?

In [119]: a=6

In [120]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [121]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [122]: class Error(Exception):
     ...:    """Base class for other exceptions"""
     ...:    error_val=4
     ...:    pass
     ...:
     ...: class ValueTooSmallError(Error):
     ...:    """Raised when the input value is too small"""
     ...:    print "LARGE ERROR"
     ...:    pass
     ...:
     ...: class ValueTooLargeError(Error):
     ...:    """Raised when the input value is too large"""
     ...:    pass
     ...:
LARGE ERROR

In [123]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [124]: class Error(Exception):
     ...:    """Base class for other exceptions"""
     ...:    error_val=4
     ...:    pass
     ...:
     ...: class ValueTooSmallError(Error):
     ...:    """Raised when the input value is too small"""
     ...:    pass
     ...:
     ...: class ValueTooLargeError(Error):
     ...:    """Raised when the input value is too large"""
     ...:    print "LARGE ERROR"
     ...:    pass
     ...:
     ...:
LARGE ERROR

In [125]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:

ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [126]: e=ValueTooLargeError()

In [127]: e
Out[127]: __main__.ValueTooLargeError()

In [128]: class Error(Exception):
     ...:    """Base class for other exceptions"""
     ...:    error_val=4
     ...:    pass
     ...:
     ...: class ValueTooSmallError(Error):
     ...:    """Raised when the input value is too small"""
     ...:    pass
     ...:
     ...: class ValueTooLargeError(Error):
     ...:    """Raised when the input value is too large"""
     ...:    return "LARGE ERROR"
     ...:
     ...:
     ...:
  File "<ipython-input-128-27ea301ca0cd>", line 12
    return "LARGE ERROR"
SyntaxError: 'return' outside function


In [129]: class Error(Exception):
     ...:    """Base class for other exceptions"""
     ...:    error_val=4
     ...:    pass
     ...:
     ...: class ValueTooSmallError(Error):
     ...:    """Raised when the input value is too small"""
     ...:    pass
     ...:
     ...: class ValueTooLargeError(Error):
     ...:    """Raised when the input value is too large"""
     ...:    def __init__(self):
     ...:        Exception.__init__(self,"LARGE ERROR")
     ...:
     ...:
     ...:
     ...:

In [130]: e=ValueTooLargeError()

In [131]: e
Out[131]: __main__.ValueTooLargeError('LARGE ERROR')

In [132]: raise e
---------------------------------------------------------------------------
ValueTooLargeError                        Traceback (most recent call last)
<ipython-input-132-0309b1d6b663> in <module>()
----> 1 raise e

ValueTooLargeError: LARGE ERROR

In [133]: e=ValueTooSmallError()

In [134]: e
Out[134]: __main__.ValueTooSmallError()

In [135]: raise e
---------------------------------------------------------------------------
ValueTooSmallError                        Traceback (most recent call last)
<ipython-input-135-0309b1d6b663> in <module>()
----> 1 raise e

ValueTooSmallError:

In [136]: e=ValueTooLargeError()

In [137]: raise e
---------------------------------------------------------------------------
ValueTooLargeError                        Traceback (most recent call last)
<ipython-input-137-0309b1d6b663> in <module>()
----> 1 raise e

ValueTooLargeError: LARGE ERROR

In [138]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
LARGE ERROR
ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?

In [139]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     raise e
     ...:     print "ERROR VALUETOOLARGE"
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
WHICH ERROR DID I PRINT?
---------------------------------------------------------------------------
ValueTooLargeError                        Traceback (most recent call last)
<ipython-input-139-e0496de048f7> in <module>()
      3         raise ValueTooLargeError
      4 except ValueTooLargeError as e:
----> 5     raise e
      6     print "ERROR VALUETOOLARGE"
      7 except IOError as e:

ValueTooLargeError: LARGE ERROR

In [140]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...: except ValueTooLargeError as e:
     ...:     print "ERROR VALUETOOLARGE"
     ...:     raise e
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
     ...:
     ...:
ERROR VALUETOOLARGE
WHICH ERROR DID I PRINT?
---------------------------------------------------------------------------
ValueTooLargeError                        Traceback (most recent call last)
<ipython-input-140-65e1037150e7> in <module>()
      4 except ValueTooLargeError as e:
      5     print "ERROR VALUETOOLARGE"
----> 6     raise e
      7 except IOError as e:
      8     print e

ValueTooLargeError: LARGE ERROR

In [141]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:         break
     ...: except ValueTooLargeError as e:
     ...:     print "ERROR VALUETOOLARGE"
     ...:     raise e
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
  File "<ipython-input-141-583c8387702b>", line 4
    break
SyntaxError: 'break' outside loop


In [142]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:     break
     ...: except ValueTooLargeError as e:
     ...:     print "ERROR VALUETOOLARGE"
     ...:     raise e
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
  File "<ipython-input-142-33b96c35fbf2>", line 4
    break
SyntaxError: 'break' outside loop


In [143]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:     break
     ...: except ValueTooLargeError as e:
     ...:     print "ERROR VALUETOOLARGE"
     ...:     raise e
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
  File "<ipython-input-143-33b96c35fbf2>", line 4
    break
SyntaxError: 'break' outside loop


In [144]: try:
     ...:     if a==6:
     ...:         raise ValueTooLargeError
     ...:     elif a==10:
     ...:         raise ValueTooSmallError
     ...:     break
     ...: except ValueTooLargeError as e:
     ...:     print "ERROR VALUETOOLARGE"
     ...:     raise e
     ...: except IOError as e:
     ...:     print e
     ...: finally:
     ...:     print "WHICH ERROR DID I PRINT?"
  File "<ipython-input-144-28012ae23aa2>", line 6
    break
SyntaxError: 'break' outside loop
__connection = connection


import logging

logger = logging.getLogger()


def f():

    try:

        flaky_func()

    except Exception:

        logger.exception()

        raise
'''
