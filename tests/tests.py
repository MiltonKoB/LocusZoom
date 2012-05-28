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

  ### 
  # TESTS THAT SHOULD FAIL:
  ###

  # Test when header is missing entirely from file. 
  header_test = Test("--metal tests/data/no_header_rs1531343.txt --refsnp rs1531343 --delim space");
  header_test.title = "Metal file missing header test";
  header_test.add_required_file("tests/data/no_header_rs1531343.txt");
  header_test.should_fail = True;
  tests.append(header_test);

  # Test when refsnp is missing from file. 
  noref_test = Test("--metal tests/data/missing_refsnp_rs1531343.txt --refsnp rs1531343 --delim space");
  noref_test.title = "Refsnp missing from metal file";
  noref_test.add_required_file("tests/data/missing_refsnp_rs1531343.txt");
  noref_test.should_fail = True;
  tests.append(noref_test);

  # Test when p-value for refsnp is missing from file. 
  nopval_test = Test("--metal tests/data/missing_refsnp_pvalue_rs1531343.txt --refsnp rs1531343 --delim space");
  nopval_test.title = "Refsnp's p-value missing from metal file";
  nopval_test.add_required_file("tests/data/missing_refsnp_pvalue_rs1531343.txt");
  nopval_test.should_fail = True;
  tests.append(nopval_test);

  # Test when delimiter is specified as something other than what is in the file. 
  baddelim_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refsnp rs8050136 --delim tab");
  baddelim_test.title = "Delimiter misspecification test";
  baddelim_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  baddelim_test.should_fail = True;
  tests.append(baddelim_test);

  # Test when data file is empty. 
  empty_test = Test("--metal tests/data/empty_file --refsnp rs1002227");
  empty_test.title = "Empty metal file test";
  empty_test.add_required_file("tests/data/empty_file");
  empty_test.should_fail = True;
  tests.append(empty_test);

  # Test when snpcol is non-standard, but --markercol specifies incorrectly. 
  snpcol_bad_test = Test("--metal tests/data/HDL_ONE_Eur_b36.tbl --refsnp rs10401969 --markercol blah --pvalcol GC.Pvalue");
  snpcol_bad_test.title = "--snpcol misspecification test";
  snpcol_bad_test.add_required_file("tests/data/HDL_ONE_Eur_b36.tbl");
  snpcol_bad_test.should_fail = True;
  tests.append(snpcol_bad_test);

  # Test when pvalcol is non-standard, but --pvalcol specifies incorrectly. 
  pvalcol_bad_test = Test("--metal tests/data/HDL_ONE_Eur_b36.tbl --refsnp rs10401969 --markercol SNPColumn --pvalcol blah");
  pvalcol_bad_test.title = "--pvalcol misspecification test";
  pvalcol_bad_test.add_required_file("tests/data/HDL_ONE_Eur_b36.tbl");
  pvalcol_bad_test.should_fail = True;
  tests.append(pvalcol_bad_test);

  # Test when hitspec is missing a SNP. 
  nosnp_hitspec_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_missing_snp.txt");
  nosnp_hitspec_test.title = "Hitspec missing SNP test";
  nosnp_hitspec_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  nosnp_hitspec_test.add_required_file("tests/data/hitlist_missing_snp.txt");
  nosnp_hitspec_test.should_fail = True;
  tests.append(nosnp_hitspec_test);

  # Test when hitspec has a malformed SNP. 
  badsnp_hitspec_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_bad_snp.txt");
  badsnp_hitspec_test.title = "Hitspec containing a bad SNP in first column";
  badsnp_hitspec_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  badsnp_hitspec_test.add_required_file("tests/data/hitlist_bad_snp.txt"); 
  badsnp_hitspec_test.should_fail = True;
  tests.append(badsnp_hitspec_test); 
  
  # Test what happens when metal file has windows line terminators. 
  winline_test = Test("--metal tests/data/windows_line_term_rs1531343.txt --refsnp rs1531343 --delim space --prefix winline");
  winline_test.title = "Windows line terminator test"
  winline_test.add_required_file("tests/data/windows_line_term_rs1531343.txt");
  tests.append(winline_test);

  # **** Testing command line argument handling. 

  # Test when chr/start are specified, but not end. 
  noend_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --chr 12 --start 64200000 --delim space");
  noend_test.title = "Specifying chr/start, but no end"
  noend_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  noend_test.should_fail = True;
  tests.append(noend_test);

  # Variation on above test - start specified, no end. 
  nostart_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --chr 12 --end 64800000 --delim space");
  nostart_test.title = "Specifying chr and end, but no start"
  nostart_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  nostart_test.should_fail = True;
  tests.append(nostart_test);

  # Failing to specify chromosome
  nochr_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --start 64200000 --end 64800000 --delim space");
  nochr_test.title = "Specifying start/end, but no chr"
  nochr_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  nochr_test.should_fail = True;
  tests.append(nochr_test);

  # Specifying flank, but flank is garbage
  badflank_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refsnp rs1531343 --flank blahhhh");
  badflank_test.title = "Flank is bad test";
  badflank_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  badflank_test.should_fail = True;
  tests.append(badflank_test);

  # User specifies an LD file to use, but file is missing. 
  missing_ld_test = Test("--metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refsnp rs1531343 --ld file_that_does_not_exist");
  missing_ld_test.title = "Missing LD file test"
  missing_ld_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  missing_ld_test.should_fail = True;
  tests.append(missing_ld_test);
  
  # User defined LD file exists, but header is missing. 
  user_ld_noheader_test = Test("--prefix user_ld_noheader_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refsnp rs7578326 --ld tests/data/bad_ld_noheader.txt --no-cleanup");
  user_ld_noheader_test.title = "User-defined LD file has no header";
  user_ld_noheader_test.should_fail = True;
  tests.append(user_ld_noheader_test);

  # User gives a refgene, but this gene does not exist in the database. 
  bad_refgene_test = Test("--prefix badrefgenetest --metal tests/data/HDL_ONE_Eur_b36_normcols.tbl --refgene BLARGHYARR --flank 100kb --no-cleanup");
  bad_refgene_test.title = "Testing when --refgene doesn't exist in database";
  bad_refgene_test.should_fail = True;
  tests.append(bad_refgene_test);

  # Invalid build is specified

  # Invalid population is specified

  # Specify flank as 0kb (should work.) 

  ###
  # TESTS THAT SHOULD GENERATE A VALID PDF: 
  ###

  # Try a metal file with p-values ranging from 0.1 to 1 (this could trip
  # an improper double -log10 transformation. 
  log_trans_test = Test("--prefix log_trans_test --metal tests/data/diagramv2_rs1531343_3MB_pval-point1to1.txt --refsnp rs1531343 --flank 100kb --no-cleanup --pvalcol P.value");
  log_trans_test.title = "Testing file with p-values [0.1,1], transform on";
  log_trans_test.should_fail = False;
  tests.append(log_trans_test);
  
  # Same metal file as before, but now with transformation turned off. 
  log_trans_off_test = Test("--prefix log_trans_off_test --no-transform --metal tests/data/diagramv2_rs1531343_3MB_pval-point1to1.txt --refsnp rs1531343 --flank 100kb --no-cleanup --pvalcol P.value");
  log_trans_off_test.title = "Testing file with p-values [0.1,1], transform off";
  log_trans_off_test.should_fail = False;
  tests.append(log_trans_off_test);
  
  # Try with a metal file that has blank lines. 
  blanklines_test = Test("--prefix blanklines_test --metal tests/data/blanklines_diagramv2_rs1531343_3MB.txt --refsnp rs1531343 --flank 300kb --delim space");
  blanklines_test.title = "Testing metal file with blank lines";
  blanklines_test.add_required_file("tests/data/blanklines_diagramv2_rs1531343_3MB.txt");
  tests.append(blanklines_test);

  # Try using a gzipped file. 
  gzip_test = Test("--prefix gzip_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt.gz --delim space --refsnp rs1531343");
  gzip_test.title = "Testing gzipped file";
  gzip_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt.gz");
  tests.append(gzip_test);
  
  # Test a metal file that has random line endings (\r\n, \r, \n.) 
  bad_lineends_test = Test("--prefix bad_lineends --metal tests/data/diagramv2_semeta_120108_bad_lineends.tbl --refsnp rs1531343 --flank 500kb --plotonly");
  bad_lineends_test.title = "Testing normal metal file with bad line endings"
  tests.append(bad_lineends_test);
  
  # Try using a gzipped with bad line endings. 
  gzip_badlineends_test = Test("--prefix gzip_badlineends --metal tests/data/diagramv2_badlineends_chr10_rs7903146_500kb.tbl --refsnp rs7903146 --flank 300kb");
  gzip_badlineends_test.title = "Testing gzipped file with bad line endings";
  gzip_badlineends_test.add_required_file("tests/data/diagramv2_badlineends_chr10_rs7903146_500kb.tbl");
  tests.append(gzip_badlineends_test);

  # Try using a bzip2 file. 
  bzip_test = Test("--prefix bzip_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt.bz2 --delim space --refsnp rs1531343");
  bzip_test.title = "Testing bz2 file";
  bzip_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt.bz2");
  tests.append(bzip_test);
  
  # Arbitrary precision test. 
  arb_prec_test = Test("--no-cleanup --prefix arb_prec_test --metal tests/data/HDL_ONE_Eur_b36.tbl --refsnp rs3764261 --markercol SNPColumn --pvalcol GC.Pvalue");
  arb_prec_test.title = "Testing arbitrary precision for extremely small p-values"
  arb_prec_test.add_required_file("tests/data/HDL_ONE_Eur_b36.tbl");
  tests.append(arb_prec_test); 

  # Find SNP automatically, with potential arb precision problems. 
  cetp_auto_test = Test("--prefix cetp_auto_test --metal tests/data/HDL_ONE_Eur_b36.tbl --refgene CETP --markercol SNPColumn --pvalcol GC.Pvalue");
  cetp_auto_test.title = "Testing auto-detection of extremely small p-value"
  cetp_auto_test.add_required_file("tests/data/HDL_ONE_Eur_b36.tbl");
  tests.append(cetp_auto_test); 

  # Try using different populations for generating LD. 
  jpt_test = Test("--prefix JPT_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop JPT+CHB");
  jpt_test.title = "Testing JPT+CHB LD info";
  tests.append(jpt_test);
  
#  hg17_jpt_test = Test("--prefix hg17_JPT_test --build hg17 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop JPT+CHB");
#  hg17_jpt_test.title = "Testing hg17 JPT+CHB LD info";
#  tests.append(hg17_jpt_test);

  yri_test = Test("--prefix YRI_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop YRI");
  yri_test.title = "Testing YRI LD info";
  tests.append(yri_test);  

#  hg17_yri_test = Test("--prefix hg17_YRI_test --build hg17 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop YRI");
#  hg17_yri_test.title = "Testing hg17 YRI LD info";
#  tests.append(hg17_yri_test);

  ceu_test = Test("--prefix CEU_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop CEU");
  ceu_test.title = "Testing CEU LD info";
  tests.append(ceu_test);
  
#  hg17_ceu_test = Test("--prefix hg17_CEU_test --build hg17 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop CEU");
#  hg17_ceu_test.title = "Testing hg17 CEU LD info";
#  tests.append(hg17_ceu_test);
  
  chrx_ceu_test = Test("--prefix CEU_test --metal tests/data/diagram+_chrx-meta_090809_nall_b36_v1.tbl --refsnp rs5945326 --pop CEU");
  chrx_ceu_test.title = "Testing CEU LD on chrX";
  tests.append(chrx_ceu_test);
  
#  chrx_hg17_ceu_test = Test("--prefix hg17_CEU_test --build hg17 --metal /home/welchr/projects/diagram+_chrx/meta/2009-09-08/diagram+_chrx-meta_090809_nall_b36_v1.tbl --refsnp rs5945326 --pop CEU");
#  chrx_hg17_ceu_test.title = "Testing hg17 CEU on chrX";
#  tests.append(chrx_hg17_ceu_test);
  
  dprime_ceu_test = Test("--prefix dprime_CEU_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refgene TCF7L2 --pop CEU ldCol=dprime");
  dprime_ceu_test.title = "Testing CEU LD with dprime";
  tests.append(dprime_ceu_test);
  
  # Test 1000G LD. 
  g1k_ceu_test = Test("--prefix 1000G_CEU_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refgene LIPC --source 1000G_Aug2009 --pop CEU");
  g1k_ceu_test.title = "Testing 1000G CEU hg18";
  tests.append(g1k_ceu_test);
  
  refsnp_g1k_ceu_test = Test("--prefix 1000G_refsnp_CEU_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refsnp chr16:55558762 --flank 300kb --source 1000G_Aug2009 --pop CEU");
  refsnp_g1k_ceu_test.title = "Testing 1000G CEU hg18, refsnp is 1000G";
  tests.append(refsnp_g1k_ceu_test);
  
  g1k_june2010_ceu_test = Test("--prefix 1000G_June2010_CEU_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refgene LIPC --source 1000G_June2010 --pop CEU");
  g1k_june2010_ceu_test.title = "Testing 1000G June2010 CEU hg18";
  tests.append(g1k_june2010_ceu_test);
  
  g1k_june2010_yri_test = Test("--prefix 1000G_June2010_YRI_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refgene LIPC --source 1000G_June2010 --pop YRI");
  g1k_june2010_yri_test.title = "Testing 1000G June2010 YRI hg18";
  tests.append(g1k_june2010_yri_test);
  
  g1k_june2010_jptchb_test = Test("--prefix 1000G_June2010_JPTCHB_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refgene LIPC --source 1000G_June2010 --pop JPT+CHB");
  g1k_june2010_jptchb_test.title = "Testing 1000G June2010 JPT+CHB hg18";
  tests.append(g1k_june2010_jptchb_test);
  
  # Testing what happens when a 1000G SNP is used with hapmap LD. This might fail..
  g1k_hapmap_test = Test("--prefix 1000G_hapmap-CEU_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refsnp chr15:56497991 --pop CEU");
  g1k_hapmap_test.title = "Testing a 1000G refSNP with hapmap LD, large region";
  tests.append(g1k_hapmap_test);
  
  # Same test as above, but smaller region that might trip a "pre-computed" lookup.
  g1k_precomp_hapmap_test = Test("--prefix 1000G_hapmap-CEU_hg18_test --metal tests/data/Lipids1kG_METAANALYSIS_HDL_FDC_GC1_hg18.tbl --refsnp chr15:56497991 --flank 100kb --pop CEU");
  g1k_precomp_hapmap_test.title = "Testing a 1000G refSNP with hapmap LD, small region";
  tests.append(g1k_precomp_hapmap_test);
  
  # Test user-defined LD file.
  user_ld_test = Test("--prefix user_ld_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --delim space --refsnp rs7578326 --ld tests/data/templd_rs7578326_and_rs1002227.txt --no-cleanup");
  user_ld_test.title = "Testing user-defined LD file";
  tests.append(user_ld_test);
  
  # Try plotting a region on X. 
  chrx_region_test = Test("--prefix chrx_region --metal tests/data/diagram+_chrx-meta_090809_nall_b36_v1.tbl --chr X --start 152000000 --end 153000000");
  chrx_region_test.title = "Testing region specifying for chr X";
  chrx_region_test.add_required_file("tests/data/diagram+_chrx-meta_090809_nall_b36_v1.tbl");
  tests.append(chrx_region_test); 

  # Try plotting a SNP on X. 
  chrx_snp_test = Test("--prefix chrx_snp --metal tests/data/diagram+_chrx-meta_090809_nall_b36_v1.tbl --refsnp rs5945326");
  chrx_snp_test.title = "Testing DIAGRAM SNP on chrX";
  tests.append(chrx_snp_test);

  # Try using a gene hitspec file. 
  t2d_genehits_test = Test("--prefix t2d_genehits --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/diagramv2_gene_hitlist.txt --delim space");
  t2d_genehits_test.title = "Testing hitspec file with only genes";
  t2d_genehits_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  t2d_genehits_test.add_required_file("tests/data/diagramv2_gene_hitlist.txt");
  tests.append(t2d_genehits_test);
  
  # Try using a hitspec that has only chr/start/stop in it. 
  t2d_chrpos_hitspec = Test("--prefix t2d_chrposhitspec --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/diagram_chrposonly_hitspec.txt --delim space");
  t2d_chrpos_hitspec.title = "Testing hitspec file with only chr/start/stop";
  t2d_chrpos_hitspec.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  t2d_chrpos_hitspec.add_required_file("tests/data/diagram_chrposonly_hitspec.txt");
  tests.append(t2d_chrpos_hitspec);
  
  # Test when hitspec has a bad line endings (either mixed up endings, multiple different endings, etc.)
  badlineends_hitspec_test = Test("--prefix badlineends_hitspec --delim space --plotonly --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_badlineends.txt");
  badlineends_hitspec_test.title = "Hitspec containing multiple line endings";
  badlineends_hitspec_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  badlineends_hitspec_test.add_required_file("tests/data/hitlist_badlineends.txt"); 
  tests.append(badlineends_hitspec_test); 
  
  # Test when hitspec has blank lines. 
  blanklines_hitspec_test = Test("--prefix blanklines_hitspec --delim space --plotonly --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_blanklines.txt");
  blanklines_hitspec_test.title = "Hitspec containing blank lines";
  blanklines_hitspec_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  blanklines_hitspec_test.add_required_file("tests/data/hitlist_blanklines.txt"); 
  tests.append(blanklines_hitspec_test); 
  
  # Every possible legal hitspec combination. 
  t2d_hitspec_allcombos_test = Test("--prefix t2d_hitspec_allcombos_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_every_combo.txt --delim space");
  t2d_hitspec_allcombos_test.title = "Testing hitspec file with every possible combo of LEGAL entries";
  t2d_hitspec_allcombos_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  t2d_hitspec_allcombos_test.add_required_file("tests/data/hitlist_every_combo.txt");
  tests.append(t2d_hitspec_allcombos_test);
  
  # Every possible illegal hitspec combination. 
  t2d_hitspec_badcombos_test = Test("--prefix t2d_hitspec_badcombos_test --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --hitspec tests/data/hitlist_every_bad_combo.txt --delim space");
  t2d_hitspec_badcombos_test.title = "Testing hitspec file with every possible combo of BAD entries";
  t2d_hitspec_badcombos_test.add_required_file("tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt");
  t2d_hitspec_badcombos_test.add_required_file("tests/data/hitlist_every_bad_combo.txt");
  tests.append(t2d_hitspec_badcombos_test);

  # Try tests where we flip SNPs from rsID to chr:pos format. 
  allchrpos_rsid = Test("--prefix allchrpos_rsid --metal tests/data/diagramv2_SEmeta_120108_chrpos_b36_v1.tbl --refsnp rs1531343 --flank 550kb --no-cleanup");
  allchrpos_rsid.title = "Region specified with rsIDs (see next plot!)";
  tests.append(allchrpos_rsid);
  
  allchrpos_chrpos = Test("--prefix allchrpos_chrpos --metal tests/data/diagramv2_SEmeta_120108_chrposONLY_b36_v1.tbl --refsnp rs1531343 --flank 550kb --no-cleanup");
  allchrpos_chrpos.title = "Region specified with only chr:pos (compare to previous)";
  tests.append(allchrpos_chrpos);

  # hg19 LD tests
  hg19_1000g_eur = Test("--prefix hg19_1000g_eur --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refgene TCF7L2 --plotonly --delim space --source 1000G_Nov2010 --pop EUR --build hg19");
  hg19_1000g_eur.title = "hg19 / 1000G_Nov2010 / EUR";
  tests.append(hg19_1000g_eur);

  hg19_1000g_afr = Test("--prefix hg19_1000g_afr --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refgene TCF7L2 --plotonly --delim space --source 1000G_Nov2010 --pop AFR --build hg19");
  hg19_1000g_afr.title = "hg19 / 1000G_Nov2010 / AFR";
  tests.append(hg19_1000g_afr);
  
  hg19_1000g_asn = Test("--prefix hg19_1000g_asn --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refgene TCF7L2 --plotonly --delim space --source 1000G_Nov2010 --pop ASN --build hg19");
  hg19_1000g_asn.title = "hg19 / 1000G_Nov2010 / ASN";
  tests.append(hg19_1000g_asn);

  ###
  # GOLD STANDARD TESTS
  # These are tests where a PDF generated is compared to a "gold standard" reference PDF, and visually compared. 
  ###

  gold1 = Test("--prefix gold1 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refsnp rs1531343 --flank 500kb --delim space --plotonly");
  gold1.title = "HMGA2 Region"
  gold1.gold_std("tests/standards/100118_rs1531343.pdf");
  tests.append(gold1);

  gold2 = Test("--prefix gold2 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refgene ZFAND6 --flank 200kb --delim space --plotonly");
  gold2.title = "ZFAND6 Region"
  gold2.gold_std("tests/standards/100118_ZFAND6.pdf");
  tests.append(gold2);

  gold3 = Test("--prefix gold3 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refgene ZFAND6 --flank 200kb --delim space --plotonly --prefix no_annot showAnnot=F");
  gold3.title = "ZFAND6 Region, without SNP annotations"
  gold3.gold_std("tests/standards/no_annot_100118_ZFAND6.pdf");
  tests.append(gold3);

  gold4 = Test("--prefix gold4 --metal tests/data/DIAGRAMv2_EU_112808_nall_results_formetal.txt --refsnp rs231362 --chr 11 --start 2348047 --end 2922200 --delim space --plotonly");
  gold4.title = "KCNQ1 Region";
  gold4.gold_std("tests/standards/100118_rs231362.pdf");
  tests.append(gold4);

  gold5 = Test("--prefix gold5 --metal tests/data/HDL_ONE_Eur_b36_normcols.tbl --refsnp rs3764261 --flank 100kb --plotonly");
  gold5.title = "CETP Super Region!!"
  gold5.gold_std("tests/standards/100118_rs3764261.pdf");
  tests.append(gold5);

  gold6 = Test("--prefix gold6 --metal tests/data/HDL_ONE_Eur_b36_normcols.tbl --refsnp rs174546 --flank 150kb --plotonly");
  gold6.title = "FADS1 Region";
  gold6.gold_std("tests/standards/100118_rs174546.pdf");
  tests.append(gold6);

  gold7 = Test("--prefix gold7 --metal tests/data/HDL_ONE_Eur_b36_normcols.tbl --refsnp rs4846914 --flank 600kb --plotonly");
  gold7.title = "GALNT2 Region"
  gold7.gold_std("tests/standards/100118_rs4846914.pdf");
  tests.append(gold7);

  gold8 = Test("--prefix gold8 --metal tests/data/HDL_ONE_Eur_b36_normcols.tbl --refsnp rs4846914 --flank 600kb --plotonly --prefix no_annot showAnnot=F");
  gold8.title = "GALNT2 Region, no SNP annotations";
  gold8.gold_std("tests/standards/no_annot_100118_rs4846914.pdf");
  tests.append(gold8);

  gold9 = Test("--metal tests/data/HDL_ONE_Eur_b36.tbl --refsnp rs10401969 --prefix gold9 --metal tests/data/HDL_ONE_Eur_b36.tbl --gene TCF7L2 --pvalcol GC.Pvalue --flank 600kb --plotonly  showAnnot=F smallDot=.5 bigDot=2 weightCol=weight");
  gold9.title = "TCF7L2 Region, dot sizes by weight";
  gold9.gold_std("tests/standards/no_annot_100118_rs4846914.pdf");
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
