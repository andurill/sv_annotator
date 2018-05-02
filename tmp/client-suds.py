import sys
from suds.client import Client

URL = 'https://mutalyzer.nl/services/?wsdl'

if len(sys.argv) < 3:
    print 'Please provide a variant'
    sys.exit(1)

c = Client(URL, cache=None)
o = c.service

print 'Checking ' + sys.argv[1] + ' ...'

r = o.numberConversion(sys.argv[1], sys.argv[2])

#if r.valid:
#    print 'Valid!'
#else:
#    print 'Not valid!'

print r

#if r.messages:
#    for m in r.messages.SoapMessage:
#        print 'Message (%s): %s' % (m.errorcode, m.message)
