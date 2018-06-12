

import six

def isIterable(x):
    return (not isinstance(x, six.string_types)) and hasattr(x, "__iter__")   