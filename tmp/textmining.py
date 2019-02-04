
import sys
from SOAPpy import WSDL

URL = 'https://mutalyzer.nl/services/?wsdl'

o = WSDL.Proxy(URL)

for line in sys.stdin.readlines() :
    for word in line.split() :
        if o.checkSyntax(variant = word).valid == 'true' :
            print word
