#!/bin/bash
 sed '23,103!d' -i  *.dat
 sed !d CJ15.txt >> *.dat
