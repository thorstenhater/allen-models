#!/bin/bash

rm -rf 473561729.zip 473561729
curl -q -o '473561729.zip' -L 'https://senselab.med.yale.edu/modeldb/eavBinDown?o=184316&a=23&mime=application/zip'
unzip 473561729.zip
