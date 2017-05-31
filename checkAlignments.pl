#!/usr/bin/env perl

use strict;
use warnings;

#Given a list of Toil alignment jobs do the following.
my $toil_cmd_file = "/gscmnt/gc2736/griffithlab_gms/Breast_cfDNA_Ademuyiwa/cwl_toil_runs/results/align_toil_cmds.sh";
my %jobs;

#1. Get the Toil log file path
#2. Get the Toil expected results dir path
#3. Get the Toil working dir path
#4. Get the Toil job store dir path

my $j = 0;
open (CMD, $toil_cmd_file) || die "\n\ncould not open Toil cmd file: $toil_cmd_file\n\n";
while(<CMD>){
  $j++;
  chomp $_;
  my $logfile = '';
  my $outdir = '';
  my $workdir = '';
  my $jobstore = '';

  if ($_ =~ /\-\-logFile\=(\S+)/){
    $logfile = $1;
  }else{
    die "\n\ncould not determine --logFile path from toil cmd:\n\n$_\n\n";
  }

  if ($_ =~ /\-\-outdir\=(\S+)/){
    $outdir = $1;
    $outdir .= "/" unless ($outdir =~ /.*\/$/);
  }else{
    die "\n\ncould not determine --outdir path from toil cmd:\n\n$_\n\n";
  }

  if ($_ =~ /\-\-workDir\=(\S+)/){
    $workdir = $1;
    $workdir .= "/" unless ($workdir =~ /.*\/$/);
  }else{
    die "\n\ncould not determine --workDir path from toil cmd:\n\n$_\n\n";
  }

  if ($_ =~ /\-\-jobStore\=(\S+)/){
    $jobstore = $1;
    $jobstore .= "/" unless ($jobstore =~ /.*\/$/);
  }else{
    die "\n\ncould not determine --jobStore path from toil cmd:\n\n$_\n\n";
  }
  $jobs{$j}{logfile} = $logfile;
  $jobs{$j}{outdir} = $outdir;
  $jobs{$j}{workdir} = $workdir;
  $jobs{$j}{jobstore} = $jobstore;
}
close(CMD);

#5. Check for success of the Toil jobs
#   - Check for both success in the Toil log and also check for existence of the results files
#If the job is successful. Do the following
#6. Delete all undesired temp files in the output directory
#7. Use `toil clean` to clean up the jobStore for this job
#8. Delete files in the working dir
#9. Write a file indicating completion of the pipeline. Then make the final results files read only.
foreach my $j (sort {$a <=> $b} keys %jobs){
  my $overall_status = "Pending";
  my $logfile = $jobs{$j}{logfile};
  my $outdir = $jobs{$j}{outdir};
  my $workdir = $jobs{$j}{workdir};
  my $jobstore = $jobs{$j}{jobstore};

  if (-e $logfile){
    my $logfile_success = 0;
    open(LOG, $logfile) || die "\n\ncould not open log file: $logfile\n\n";
    while(<LOG>){
      chomp $_;
      if ($_ =~ /Finished\s+toil\s+run\s+successfully/){
        $logfile_success = 1;
      }
    }
    close(LOG);

    my $file_success = 0;
    if (-e $outdir."final.crai" && -e $outdir."final.cram" && -e $outdir."final.cram.crai"){
      $file_success = 1;
    }

    $overall_status = "Complete" if ($logfile_success == 1 && $file_success == 1);
  }

  if ($overall_status eq "Complete"){
    my $files = $outdir . "lsf* " . $outdir . "tmp*";
    my $rm_cmd = "rm -fr $files";
    print "\n$rm_cmd";
    system($rm_cmd);

    my $toil_clean_cmd = "toil clean $jobstore";
    print "\n$toil_clean_cmd";
    system($toil_clean_cmd);

    my $workdir_rm_cmd = "rm -fr $workdir";
    print "\n$workdir_rm_cmd";
    system($workdir_rm_cmd);

    my $complete_file = $outdir . "COMPLETE";
    unless (-e $complete_file){
      my $complete_file_cmd = "touch $complete_file";
      print "\n$complete_file_cmd";
      system($complete_file_cmd);

      my $chmod_cmd = "chmod -w $outdir";
      print "\n$chmod_cmd";
      system($chmod_cmd);
    }

  }
  print "\n";
  print "$j\t$jobs{$j}{workdir}\t$overall_status\n";


}


print "\n\n";

exit;


