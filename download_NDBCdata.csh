#!/bin/csh 

set year = 2017
set d0 = 'latest'
set d1 = '2017-01-01T00:00Z'
set d2 = '2017-01-15T23:59Z'
set d  = "${d1}/${d2}"
rm wavesfile.csv
set url = 'https://sdf.ndbc.noaa.gov/sos/server.php'
set request = 'request=GetObservation&service=SOS&version=1.0.0&offering=urn:ioos:network:noaa.nws.ndbc:all&observedproperty=waves&responseFormat=text/csv'
wget -c -r -l 1 -np "${url}?${request}&eventtime=${d}" -O wavesfile.csv

set filename = 'wavesfile.csv'



# call python script to read data into an array
python readNDBC.py  $filename >>& NDBC_read.log 
