#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2010 Ryan Welch, Randall Pruim
#
# LocusZoom is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LocusZoom is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import os, sys, re, gzip, urllib, urllib2, tempfile, ftplib, time
import os.path as path
from optparse import OptionParser

# Constants.
DBMEISTER = "bin/dbmeister.py";
SQLITE_SNP_POS = "snp_pos";
SQLITE_TRANS = "refsnp_trans";
SQLITE_REFFLAT = "refFlat";
SQLITE_GENCODE = "gencode";
RS_MERGE_ARCH_URL = "ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606/RsMergeArch.bcp.gz";
GENCODE_URL = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/";
GENCODE_FTP = "ftp.sanger.ac.uk";
GUNZIP_PATH = "gunzip";
REFFLAT_HEADER = "geneName name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds".split();

UCSC_TO_GRC = {
  'hg19' : 'GRCh37'
}

# Utility functions.

def which(f):
  for path in os.environ['PATH'].split(os.pathsep):
    if not os.path.exists(path):
      continue;

    if not os.path.isdir(path):
      continue;

    for file in os.listdir(path):
      if os.path.basename(f) == file:
        return os.path.join(path,file);

  return None;

def find_relative(file):
  full_path = None;

  # Find the m2zfast root, using the script's location.
  start_loc = os.path.realpath(sys.argv[0]);
  script_dir = None;
  if os.path.isdir(start_loc):
    script_dir = start_loc;
  else:
    script_dir = os.path.dirname(start_loc);
  root_dir = os.path.join(script_dir,"../");

  # If the file to find has a path, it means it is a path relative
  # to the m2zfast root. We need to attach that path to the root.
  (file_path,file_name) = os.path.split(file);
  if file_path != "":
    root_dir = os.path.abspath(os.path.join(root_dir,file_path));

  if file_name == "" or file_name is None:
    if os.path.exists(root_dir):
      full_path = root_dir;
  else:
    temp_path = os.path.join(root_dir,file_name);
    if os.path.exists(temp_path):
      full_path = temp_path;

  return full_path;

# Tries to find a file in the following order:
# 1) the given file (duh)
# 2) relative to m2zfast's root directory
# 3) on the user's path
def find_systematic(file):
  if file is None:
    return None;

  if os.path.isfile(file):
    return os.path.abspath(file);

  relative = find_relative(file);
  if relative:
    return relative;

  whiched_file = which(file);
  if whiched_file is not None:
    return whiched_file;

  return None;

# Function to see if a given URL actually exists.
def exists(url):
  try:
    urllib2.urlopen(url);
  except urllib2.HTTPError, e:
    if e.code == 404:
      return False;

  return True;

def dl_hook(count,block_size,total_size):
  percent = min(100,count*block_size*100.0/total_size);
  sys.stdout.write("\r%.1f%%" % percent)
  sys.stdout.flush()

class UCSCManager:
  UCSC_MAIN_URL = "http://hgdownload.cse.ucsc.edu/goldenPath";
  UCSC_FTP_URL = "ftp://hgdownload.cse.ucsc.edu/goldenPath/";

  def __init__(self):
    pass

  @staticmethod
  def getLatestHumanBuild():
    latest_hg = None;

    resp = urllib2.urlopen(UCSCManager.UCSC_FTP_URL);
    lines = resp.readlines();
    dirs = [i.rstrip().split()[8] for i in lines];
    p = re.compile("hg(\d+)");
    hg = filter(lambda x: p.search(x) is not None,dirs);
    hg_versions = map(lambda x: int(p.search(x).groups()[0]),hg);
    latest_hg = sorted(hg_versions,reverse=True)[0];

    return "hg" + str(latest_hg);

  @staticmethod
  def getLatestSNPTable(build):
    p = re.compile("snp(\d+?).sql");
    resp = urllib.urlopen(UCSCManager.UCSC_MAIN_URL + "/" + build + "/" + "database");
    tables = set();
    for line in resp:
      m = p.search(line);
      if m is not None:
        table = "snp" + str(m.groups()[0]);
        tables.add(table);

    return max(tables);

  @staticmethod
  def downloadLatestSNPTable(dir,build):
    latest_table = UCSCManager.getLatestSNPTable(build);

    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database',latest_table + '.txt.gz']);

    file = path.join(dir,latest_table + ".gz");
    #progress = urlgrabber.progress.TextMeter();
    #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);

    urllib.urlretrieve(url,file);
    #grabber.urlgrab(url,file);

    return file;

  @staticmethod
  def downloadLatestRefFlat(dir,build):
    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database','refFlat.txt.gz']);
    file = path.join(dir,'refFlat_' + build + '.txt.gz');

    #progress = urlgrabber.progress.TextMeter();
    #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);
    urllib.urlretrieve(url,file,reporthook=dl_hook);
    #grabber.urlgrab(url,file);

    return file;

  @staticmethod
  def download_snp_table(dir,build,table):
    url = "/".join([UCSCManager.UCSC_MAIN_URL,build,'database',table + '.txt.gz']);
    file = path.join(dir,table + ".gz");

    if not exists(url):
      print >> sys.stderr, "Could not find SNP table %s at UCSC - check your table name." % table;
      print >> sys.stderr, "URL attempted was: %s" % url;
      sys.exit(1);

    try:
      urllib.urlretrieve(url,file,reporthook=dl_hook);
      #progress = urlgrabber.progress.TextMeter();
      #grabber = urlgrabber.grabber.URLGrabber(progress_obj=progress,timeout=30);
      #grabber.urlgrab(url,file);
    except IOError:
      print >> sys.stderr, "A network connection to the UCSC data repository could not be made.";
      sys.exit(1);

    return file;

class MergeHistory:
  def __init__(self):
    self.graph = {};

  def add_merge(self,source,target):
    self.graph[source] = target;

  def get_merge_target(self,source):
    return self.graph.get(source);

  def find_current(self,source):
    target = source;
    while 1:
      next_t = self.get_merge_target(target);
      if next_t is not None:
        target = next_t;
      else:
        break;

    return target;

  def iter_node(self):
    for node in self.graph:
      yield node;

def parse_rsmerge(file,snp_set,snp_build):
  snp_build = int(snp_build.replace("snp",""));
  if os.path.splitext(file)[1] == ".gz":
    f = gzip.open(file);
  else:
    f = open(file);

  hist = MergeHistory();
  for line in f:
    e = line.rstrip().split("\t");
    build = int(e[2]);

    if build > snp_build:
      continue;

    rs_high = "rs" + e[0];
    rs_low = "rs" + e[1];

    hist.add_merge(rs_high,rs_low);

  refsnp_trans = {};
#  for snp in hist.iter_node():
#    if snp not in snp_set:
#      latest = hist.find_current(snp);
#      refsnp_trans[snp] = latest;

  for snp in hist.iter_node():
    latest = hist.find_current(snp);
    refsnp_trans[snp] = latest;

  return refsnp_trans;

def download_merge_arch(merge_arch_file="RsMergeArch.bcp.gz"):
  if os.path.isfile(merge_arch_file):
    os.remove(merge_arch_file);

  urllib.urlretrieve(RS_MERGE_ARCH_URL,merge_arch_file,reporthook=dl_hook);
  print "";

  return merge_arch_file;

def download_gencode(gencode_release,build):
  """
  Download the latest annotation GTF from GENCODE for a particular release.
  Also makes an attempt to check this directory for the correct genome build.
  """

  basedir = "pub/gencode/Gencode_human/release_%s" % gencode_release;

  # Check that the directory has files corresponding to the correct genome build.
  ftp = ftplib.FTP(GENCODE_FTP);
  ftp.login("anonymous","locuszoom");
  ftp.cwd(basedir);

  files = [];
  try:
    files = ftp.nlst()
  except ftplib.error_perm, resp:
    if str(resp) == "550 No files found":
      pass
    else:
      print >> sys.stderr, "Error: could not list files in GENCODE FTP, skipping GENCODE...";
      return None;

  # Check for the correct genome build.
  grc_build = UCSC_TO_GRC.get(build);
  if grc_build is not None:
    # Do we have the GRC build in the genome.fa file?
    build_ok = any(map(lambda x: re.search(r'.*%s.*\.fa\.gz' % grc_build,x) is not None,files));

    # If we couldn't find a GRCh## FA file that matches the correct build, we shouldn't proceed.
    if not build_ok:
      print >> sys.stderr, "Error: GENCODE release %s appears to not match the genome build (%s/%s) that you requested." % (gencode_release,build,grc_build);
      return None;
  else:
    print >> sys.stderr, "Warning: could not convert build %s to a GRC build (e.g. GRCh37). This means skipping the" \
                         " check that your UCSC genome build is valid for this GENCODE release. Be absolutely sure that" \
                         " the version you're requesting matches the correct genome build." % build;

  # Now that we've done the genome build check, download the annotation file.
  url = GENCODE_URL + "release_{release}/gencode.v{release}.annotation.gtf.gz".format(release = gencode_release);
  dlfile = "gencode.v{release}.annotation.gtf.gz".format(release = gencode_release);
  urllib.urlretrieve(url,dlfile,reporthook=dl_hook);

  return dlfile;

def load_snp_set(snp_pos_file):
  snp_set = set();
  with open(snp_pos_file) as f:
    header = f.readline().split();
    snp_col = header.index("snp");

    for line in f:
      snp = line.split()[snp_col];
      snp_set.add(snp);

  return snp_set;

def fix_refflat(filepath,outpath):
  """
  Performs the following actions to a downloaded refFlat.gz file:
  1) Add header row
  2) Remove alternative haplotype and other random chromosomes
  3) Write out to new tab-delimited file
  """

  if filepath.endswith(".gz"):
    f = gzip.open(filepath);
  else:
    f = open(filepath);

  outpath = outpath.replace(".gz","");
  with f, open(outpath,'w') as out:
    print >> out, "\t".join(REFFLAT_HEADER);
    for line in f:
      e = line.split("\t");
      if "_" in e[2]:
        continue;

      out.write(line);

  return outpath;

def write_refsnp_trans(rs_merge_file,snp_set,snp_table_version,out_file):
  snp_trans = parse_rsmerge(rs_merge_file,snp_set,snp_table_version);

  # add in SNPs that didn't change names
  for snp in snp_set:
    if snp not in snp_trans:
      snp_trans[snp] = snp;

  with open(out_file,"w") as out:
    print >> out, "rs_orig\trs_current";

    for orig,cur in snp_trans.iteritems():
      print >> out, "%s\t%s" % (orig,cur);

class ChromConverter:
  def __init__(self):
    self.pattern = re.compile("(chrom|chr|)(\w+)");

  def __call__(self,chr):
    if chr is None:
      return None;

    chr = str(chr);
    search = self.pattern.search(chr);
    chr_int = None;
    if search is not None:
      (chr_string,chr_val) = search.groups();
      if chr_val is not None:
        try:
          chr_int = int(chr_val);
        except:
          if chr_val == 'X':
            chr_int = 23;
          elif chr_val == 'Y':
            chr_int = 24;
          elif chr_val == 'mito':
            chr_int = 25;
          elif chr_val == 'XY':
            chr_int = 26;
          else:
            chr_int = None;

    return chr_int;

chrom2chr = ChromConverter();

def parse_ucsc_snp_table(filepath,out):
  SNP_COL = 4;
  CHROM_COL = 1;
  POS_COL = 3;

  if filepath[-3:] == ".gz":
    f = gzip.open(filepath);
  else:
    f = open(filepath);

  out = open(out,"w");

  snp_set = set();
  with out, f:
    print >> out, "snp\tchr\tpos";

    for line in f:
      e = line.strip().split("\t");
      (snp,chr,pos) = [e[i] for i in (SNP_COL,CHROM_COL,POS_COL)];

      snp_set.add(snp);

      # fix chrom
      chr_fixed = chrom2chr(chr);
      if chr_fixed is None:
        #print >> sys.stderr, "Skipping %s, chrom invalid: %s" % (str(snp),str(chr));
        continue;

      print >> out, "\t".join(map(str,[snp,chr_fixed,pos]));

  return snp_set;

def get_settings():
  p = OptionParser();
  p.add_option("-b","--build",help="Genome build (UCSC convention), e.g. hg18, hg19, etc.",default="hg19");
  p.add_option("--gencode",help="Also build a gene table using GENCODE. This specifies the relase number.");
  p.add_option("--db",help="Database name. Defaults to locuszoom_%build%.db.")
  p.add_option("--no-cleanup",help="Leave temporary files alone after creating database instead of deleting them.",default=False,action="store_true");

  opts, args = p.parse_args();

  # if opts.tmpdir is None:
  #   opts.tmpdir = tempfile.mkdtemp();

  if opts.db is None:
    opts.db = "locuszoom_%s.db" % opts.build;

  if opts.gencode is not None:
    # Get rid of the "v" if they accidentally specify it.
    opts.gencode = opts.gencode.replace("v","");

    # Should be an integer.
    try:
      int(opts.gencode);
    except:
      raise ValueError, "Error: --gencode should specify a release number (e.g. 17, 18, 19, ...)";

  return opts, args;

def mkpath(tmpdir,filename):
  return os.path.join(tmpdir,filename);

def main():
  opts, args = get_settings();
  genome_path = ".";

  # If no build was specified, find the latest human genome build from UCSC.
  build = opts.build;
  if build is None:
    build = UCSCManager.getLatestHumanBuild();

  # If we already have the latest SNP table for this build, do nothing.
  # Otherwise, grab the latest table from UCSC and install it.
  print "Asking UCSC for latest SNP table in build %s.." % build;
  ucsc_latest_snp_table = UCSCManager.getLatestSNPTable(build);
  print "Found: %s" % ucsc_latest_snp_table;

  print "Downloading SNP table %s.." % ucsc_latest_snp_table;
  snp_table_file = UCSCManager.download_snp_table(genome_path,build,ucsc_latest_snp_table);
  print "\nFinished downloading %s.." % ucsc_latest_snp_table;

  # Parse out SNP table into needed columns.
  print "Reformatting SNP table for insertion into database..";
  fixed_snp_tab = "snp_pos_%s.tab" % ucsc_latest_snp_table;
  snp_set = parse_ucsc_snp_table(snp_table_file,fixed_snp_tab);

  print "Downloading refFlat table..";
  refflat = UCSCManager.downloadLatestRefFlat(genome_path,build);
  print "\nFinished downloading %s.." % refflat;

  print "Reformatting refFlat..";
  fixed_refflat = fix_refflat(refflat,refflat.replace(".txt.gz",".tab"));

  print "Downloading RsMergeArch..";
  merge_file = download_merge_arch();

  print "Creating refsnp_trans file..";
  refsnp_trans_file = "refsnp_trans.tab";
  write_refsnp_trans(merge_file,snp_set,ucsc_latest_snp_table,refsnp_trans_file);

  if opts.gencode is not None:
    print "Downloading GENCODE file..";
    gencode_file = download_gencode(opts.gencode,build);
    print "\nFinished downloading %s.." % gencode_file;

  # Insert the files created above into the database.
  print "Creating database: %s" % opts.db;
  db_script = find_systematic(DBMEISTER);
  db_name = opts.db;
  os.system("%s --db %s --snp_pos %s" % (db_script,db_name,fixed_snp_tab));
  os.system("%s --db %s --refflat %s" % (db_script,db_name,fixed_refflat));
  os.system("%s --db %s --trans %s" % (db_script,db_name,refsnp_trans_file));

  if opts.gencode is not None:
    os.system("%s --db %s --gencode %s" % (db_script,db_name,gencode_file));

  # Do we have recombination rates for this build?
  recomb_file = find_systematic("data/build/%s/recomb_rate/recomb_rate.tab" % build);
  if recomb_file is not None:
    os.system("%s --db %s --recomb_rate %s" % (db_script,db_name,recomb_file));
  else:
    print >> sys.stderr, "Could not find a recombination rate file for this genome build.";

  # Do we have a SNP set file for this build?
  snp_set_file = find_systematic("data/build/%s/snp_set/snp_set.tab" % build);
  if snp_set_file is not None:
    os.system("%s --db %s --snp_set %s" % (db_script,db_name,snp_set_file));

  # Write a file so we know when the database was created.
  db_info = opts.db + ".info";
  with open(db_info,'w') as info_out:
    print >> info_out, time.strftime("Database created at %H:%M:%S %Z on %B %d %Y");

  # Delete all of the temporary files/directories we created.
  if not opts.no_cleanup:
    for f in [snp_table_file,fixed_snp_tab,refflat,fixed_refflat,merge_file,refsnp_trans_file,gencode_file]:
      try:
        os.remove(f);
      except:
        pass

  print "\nDatabase successfully created: %s" % db_name;
  print "To use this database, you can either: \n" \
        "1) pass it to locuszoom using the --db argument, \n" \
        "2) overwrite an existing database in <lzroot>/data/database/ (backup the existing one first!), or \n" \
        "3) add it wherever you would like, and then modify <lzroot>/conf/m2zfast.conf to point to it\n";

if __name__ == "__main__":
  main();

