#!/usr/bin/env python
import sys

#===============================================================================
# DEPRECATED / NO LONGER IN USE
# Python 2.6.x bugfix release fixed issue with gzip.open(..,"rU")
#===============================================================================

try:
  import gzip
  
  def gz_univ_readline(file):
    f = gzip.open(file);
    size = 1024 * 1024;
    read = f.read;
    linesep = "\n";
  
    buf = read(size);
    while buf: 
      buf = buf.replace("\r\n",linesep);
      buf = buf.replace("\r",linesep);
      buf = buf.replace("\n",linesep);
  
      for line in buf.split("\n"):
        yield line;
    
      buf = read(size);
  
    f.close(); 
except NameError: 
  def gz_univ_readline(file):
    msg = "Error: gzip is not available on your system. Please extract your "\
          "--metal file first before passing into LocusZoom.";
    print >> sys.stderr, msg;
    sys.exit(1);
