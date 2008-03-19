#!/usr/bin/python
import sys
from lxml import etree
from copy import deepcopy

for action, element in etree.iterparse(sys.stdin, tag="dm"):
  dm_copy=deepcopy(element)
  print dm_copy.xpath('count(row/entry[ number(.) < 0.1 ])')

