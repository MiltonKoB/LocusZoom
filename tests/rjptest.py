
#!/usr/bin/env python
import os
import sys
import re
import time
import hashlib
import tempfile
from killableprocess import *

PYEXE = "python";

# Fix path of script to be absolute.
sys.argv[0] = os.path.abspath(sys.argv[0]);

# Import m2z modules
sys.path.append(
  os.path.join(os.path.dirname(sys.argv[0]),"../src")
)

from m2zutils import *

# Import testing modules
sys.path.append(
  os.path.join(os.path.dirname(sys.argv[0]),"src/")
)

from pyPdf import PdfFileWriter, PdfFileReader 

# Force data to be written to file. 
def hard_write(f,data):
  print >> f, data;
  f.flush();
  os.fsync(f);

# Combine a list of PDFs into a single file. 
# Only combines the first page of each PDF. 
def combine_plots(pdf_list,out_file):
  output = PdfFileWriter();
  for pdf in pdf_list:
    if not os.path.isfile(pdf):
      continue;
  
    pdfobj = PdfFileReader(file(pdf,"rb"));
    output.addPage(pdfobj.getPage(0));

  out_stream = file(out_file,"wb");
  output.write(out_stream);
  out_stream.close();

# Given a directory, give a list of all files 
# in all directories underneath. 
def walk_flat(path):
  files = [];
  for i in os.walk(path):
    real_path = os.path.realpath(i[0]);
    for f in i[2]:
      full_file = os.path.join(real_path,f);
      files.append(full_file);
  
  return files;

# Get the size of a file as a long. 
def getFileSize(file):
  return os.stat(file).st_size;

# Filter a list of files down to only those with a PDF file extension. 
def filterPDF(file_list):
  p = re.compile("\.pdf",re.I);
  return filter(lambda x: p.search(x) != None,file_list);

# Class for watching a directory tree for file creations/deletions/changes. 
class FileSystemWatcher():
  def __init__(self,dir,recursive=False):
    self.dir = dir;
    self.snap_start = None;
    self.snap_end = None;
    self.recursive = recursive;
  
    if not os.path.isdir(dir):
      raise ValueError, "Error: dir %s does not exist." % str(dir);

  # Start watching the file system below self.dir for changes. 
  def start(self):
    if self.recursive:
      self.snap_start = walk_flat(self.dir);
    else:
      self.snap_start = os.listdir(self.dir);
 
    self.snap_end = None; # clear out end snapshot
  
  # Stop watching. 
  def end(self):
    if self.recursive:
      self.snap_end = walk_flat(self.dir);
    else:
      self.snap_end = os.listdir(self.dir);

  # Find files that were created between start/end. 
  def find_new(self):
    return sorted(list(set(self.snap_end).difference(set(self.snap_start))));
      
  # Find files that were deleted between start/end. 
  def find_removed(self):
    return sorted(list(set(self.snap_start).difference(set(self.snap_end))));

  # Find files whose modification times changed between start/end. 
  def find_mod_changed(self):
    pass

# Class responsible for running a set of Tests. 
class TestSuite():
  def __init__(self,bin,args,run_dir,log_file=None):
    self.bin = find_systematic(bin);
    
    if args != None:
      self.bin += " " + " ".join(args);
    
    self.run_dir = find_systematic(run_dir);
    self.tests = [];
    self.log = log_file;
    if not os.path.isdir(self.run_dir):
      raise ValueError, "Error: dir %s is not a directory." % str(run_dir);

  def add_test(self,test):
    test.id = len(self.tests);
    self.tests.append(test);

  def run_all(self):
    cwd = os.getcwd();
    
    # Try to change to run directory. 
    os.chdir(self.run_dir);
 
    # Create directory to work within.
    work_dir = time.strftime("%Y-%m-%d_%H-%M-%S");
    os.mkdir(work_dir);
    os.chdir(work_dir); 

    # Setup log file if requested. 
    log_file = None;
    if self.log:
      log_file = open(self.log,"a");
    else:
      log_file = sys.stdout;

    # Run tests. 
    test_results = {};
    all_pdfs = [];
    for test in self.tests:
      hard_write(log_file,"$ Executing test [%i]: %s" % (test.id,test.title));
      
      if test.bin == None:
        test.bin = self.bin;
     
      if log_file != None:
        test.out_file = log_file;
      
      result = test.run();

      for p in test.pdfs:
        all_pdfs.append(p);        

      for r in result:
        test_results.setdefault(test.id,[]).append(r);

    # Combine PDFs. 
    out_name = "allpdfs_" + time.strftime("%Y-%m-%d_%H-%M-%S") + ".pdf";
    combine_plots(all_pdfs,out_name);

    # Combine gold standard plots. 
    gold_list = [];
    for test in self.tests:
      if test.gold_standard != None and test.gold_standard != "":
        gold_list.append(test.gold_standard);
        for p in test.pdfs:
          gold_list.append(p);

    if len(gold_list) > 0:
      gold_out = "gold_" + time.strftime("%Y-%m-%d_%H-%M-%S") + ".pdf";
      combine_plots(gold_list,gold_out);

    # Change back to original directory. 
    os.chdir(cwd);

    # Write out summary of test results. 
    for i in xrange(len(test_results)):
      result = test_results[i];
      test = self.tests[i];

      hard_write(log_file,"Test [%i] - %s:" % (i,test.title));
      for r in result:
        hard_write(log_file,"Result:");
        hard_write(log_file,"|-  Pass: [%s]" % str(r['pass']));
        hard_write(log_file,"|-  Message: %s" % str(r['msg']));
        hard_write(log_file,"|-  File: %s" % ("NA" if not r.has_key('file') else str(r['file'])));
        hard_write(log_file,"-------");
  
      hard_write(log_file,"");

class Test():
  def __init__(self,cmd_string,out_file=None,timeout=-1):
    # Set before execution:
    self.id = None; 
    self.title = "";
    self.bin = None;
    self.cmd = cmd_string;
    self.required_files = [];
    self.timeout = timeout;
    self.gold_standard = "";
    self.should_fail = False;
    self.out_file = out_file;

    # Set after execution: 
    self.files = [];
    self.pdfs = [];

  def add_required_file(self,file):
    self.required_files.append(find_systematic(file));

  def gold_std(self,file):
    self.gold_standard = find_systematic(file);

  def hash(self):
    return hashlib.md5(self.cmd + self.title + "".join(required_files));

  def run(self):
    self.files = [];
    self.pdfs = [];

    # Make sure required files exist. 
    for f in self.required_files:
      if not os.path.isfile(f):
        return [{'pass':False,'msg':"Required file missing: %s" % str(f)}];
    
    # Form command line. 
    run_cmd = None;
    if self.bin:
      run_cmd = PYEXE + " " + self.bin + " " + self.cmd;
    else:
      raise ValueError, "Bin not specified.";

    # Should we be writing output to a file? 
    out = sys.stdout;
    if self.out_file != None:
      if isinstance(self.out_file,str):
        out = open(self.out_file,"a");
      elif isinstance(self.out_file,file):
        out = self.out_file;
      else:
        raise ValueError, "Test out_file parameter must either be a filename (string) or file object.";

    # Was there a title? 
    if self.title != "":
      run_cmd += " title=\"%s\"" % str(self.title);
      
    print run_cmd;

    # Start watching directories below the run directory for changes. 
    watcher = FileSystemWatcher(dir=".",recursive=True);
    watcher.start();

    # Run the test. 
    proc_stream = tempfile.NamedTemporaryFile("r+");
    proc = Popen(run_cmd,shell=True,stdout=proc_stream,stderr=proc_stream);
    retcode = proc.wait(timeout=self.timeout);
    watcher.end();

    # Get output from process. 
    proc_stream.seek(0);
    proc_string = proc_stream.read();
    proc_stream.close();    

    # Test results stored here. 
    results = []; 

    # Did the program catch an error? 
    caught_error = re.compile("error",re.I).search(proc_string);
    if caught_error != None:
      if self.should_fail:
        results.append({'pass':True,'msg':"Error was detected by program."});
      else:
        results.append({'pass':False,'msg':"An error was detected by the program."});
    else:
      if self.should_fail:
        results.append({'pass':False,'msg':"Test was expected to fail, but no error message was detected."});
        
    # Were any tracebacks generated? 
    traceback = re.compile("Traceback").search(proc_string);
    if traceback != None:
      if self.should_fail:
        results.append({'pass':False,'msg':"Test failed as expected, but traceback was detected."});
      else:
        results.append({'pass':False,'msg':"A traceback was detected."});        

    # Check for PDF creation. 
    self.files = watcher.find_new();
    self.pdfs = filterPDF(self.files);

    if len(self.pdfs) == 0:
      if not self.should_fail:
        results.append({'pass':False,'msg':"No PDF generated."});
    else:
      for pdf in self.pdfs:
        size = getFileSize(self.pdfs[0]);
        if size > 0:
          results.append({'pass':True,'msg':"PDF success.",'file':pdf});
        else:
          results.append({'pass':False,'msg':"PDF generated, but file size was 0.",'file':file});

    # Redirect output. 
    if self.out_file != None:
      hard_write(out,proc_string);
    else:
      out.write(proc_string);
 
    return results;

def create_tests():
  print "Creating test cases..";
  
  tests = []; 

  gold9 = Test("--metal tests/data/HDL_ONE_Eur_b36.tbl --prefix gold9 --metal tests/data/HDL_ONE_Eur_b36.tbl --refgene TCF7L2 --markercol SNPColumn --pvalcol GC.Pvalue --flank 60kb showAnnot=FALSE smallDot=0.8 largeDot=2.0 weightCol=Weight weightRange=90000,100000");
  gold9.title = "TCF7L2 Region, dot sizes by weight";
  gold9.gold_std("tests/standards/chr10_114639998-114977426.pdf");
  tests.append(gold9);
  return tests;

def main():
  test_list = create_tests();
  
  print "Setting up tests..";
  suite = TestSuite(
    bin="src/m2zfast.py",
    args=[],
    run_dir="tests/results",
    log_file="log_file.txt"
   );

  print "Adding test cases to testing suites..";
  for test in test_list:
    suite.add_test(test);

  print "Executing tests..";
  suite.run_all();

if __name__ == "__main__":
  main();
