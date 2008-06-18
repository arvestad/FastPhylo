#!/usr/bin/python
import sys
from lxml import etree
from copy import deepcopy

maxcount=0
for action, element in etree.iterparse(sys.stdin, tag="run"):
  run_copy=deepcopy(element)
  count=int(dm_copy.xpath('count'))
  if ( count > maxcount ):
    maxcount=count
print maxcount
