#!/usr/bin/env genome-perl

use warnings;
use strict;
use Data::Dumper;

#NOTE: This script will create two bash files of toil commands to be run for the breast cancer cf-dna project
#To run these toil commands you will have to set up the correct environment to launch them.  Two steps:
#STEP 1. 
#Get a *clean* environment using a docker-interactive session along with a Docker Image that contains *Toil*, ex. mgibio/unaligned-to-bqsr
#This docker image should also be the one that contains tools needed for the pipeline you want to run 
#e.g. for alignment (mgibio/unaligned-to-bqsr), or somatic variant calling (mgibio/cle), or rna-seq (mgibio/rnaseq)

#LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q docker-interactive -a 'docker(mgibio/unaligned-to-bqsr)' -Is $SHELL
#export LSB_SUB_ADDITIONAL="docker(mgibio/unaligned-to-bqsr)"

#LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q docker-interactive -a 'docker(mgibio/cle)' -Is $SHELL
#export LSB_SUB_ADDITIONAL="docker(mgibio/cle)"

#STEP 2. Set environment variables needed within that docker session:
#export LSB_DEFAULTQUEUE="research-hpc"
#export LSF_SERVERDIR="/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc"
#export LSF_LIBDIR="/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/lib"
#export LSF_BINDIR="/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin"
#export LSF_INCLUDEDIR="/opt/lsf9/9.1/include"
#export LSF_ENVDIR="/opt/lsf9/conf"
#export PATH="/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin:$PATH"

#For convenience these will be written to accompanying BASH files

#INPUTS
my %assets;
my $base_dir = "/gscmnt/gc2736/griffithlab_gms/Breast_cfDNA_Ademuyiwa/cwl_toil_runs"; $assets{$base_dir} = "d";
my $working_dir = "/gscmnt/gc2736/griffithlab_gms/Breast_cfDNA_Ademuyiwa/cwl_toil_runs/results/"; $assets{$working_dir} = "d";
my $alignment_dir = $working_dir . "alignment/"; $assets{$alignment_dir} = "d";
my $somatic_dir = $working_dir . "somatic/"; $assets{$somatic_dir} = "d";
my $sample_data_file = $working_dir . "All-Data-Summary-v5.tsv"; $assets{$sample_data_file} = "f";
my $lims_data_file = $working_dir . "LIMS_InstrumentDataInfo.tsv"; $assets{$lims_data_file} = "f";
my $sample_legend_file = $working_dir . "Sample-Legend.tsv"; $assets{$sample_legend_file} = "f";
my $alignment_yml_name = "unaligned_bam_to_bqsr.yml";

my $reference = "/gscmnt/gc2764/cad/HCC1395/arvados/refseq/GRCh38DH/GRCh38_full_analysis_set_plus_decoy_hla.fa"; $assets{$reference} = "f";
my $dbsnp = "/gscmnt/gc2764/cad/HCC1395/arvados/refseq/GRCh38DH/Homo_sapiens_assembly38.dbsnp138.vcf.gz"; $assets{$dbsnp} = "f";
my $known_indels = "/gscmnt/gc2764/cad/HCC1395/arvados/refseq/GRCh38DH/Homo_sapiens_assembly38.known_indels.vcf.gz"; $assets{$known_indels} = "f";
my $mills = "/gscmnt/gc2764/cad/HCC1395/arvados/refseq/GRCh38DH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; $assets{$mills} = "f";

my $align_cwl_exome = $base_dir . "/cancer-genomics-workflow/unaligned_bam_to_bqsr/workflow.cwl"; $assets{$align_cwl_exome} = "f";
my $align_cwl_swift = $base_dir . "/cancer-genomics-workflow/unaligned_bam_to_bqsr/workflow_no_dup_marking.cwl"; $assets{$align_cwl_swift} = "f";
$assets{"$base_dir/workDir"} = "d";
$assets{"$base_dir/jobStore"} = "d";

#Use a job group to controls slots used by these jobs
#To create a job group with 5 job limit: bgadd -L 5 /mgtoil
#To check this job group: bjobs -g /mgtoil
#To see details on this job group: bjgroup -s /mgtoil
#To adjust job limit: bgmod -L 10 /mgtoil

my $job_group = "/mgtoil";

#OUTPUTS
my $align_bash_file = $working_dir . "align_toil_cmds.sh";
my $align_session_file = $working_dir . "align_toil_session.txt";
my $somatic_bash_file = $working_dir . "somatic_toil_cmds.sh";
my $somatic_session_file = $working_dir . "somatic_toil_session.txt";

#Load data file with one line per piece of instrument data for this project
my $sample_data = &loadData('-data_file'=>$sample_data_file, '-assets'=>\%assets);

#Determine tumor/normal pairs for comparison
my $pair_data = &determinePairs('-sample_data'=>$sample_data);

print Dumper $pair_data;

#Get library and flowcell info for each piece of instrument data
#Need to create something like this for each instrument data input:
# '@RG\tID:2895804954\tPL:ILLUMINA\tPU:H75J7BBXX-CATTTTAT-CATTTTAT.8\tLB:H_UJ-GPD1001_NTN001-M1502596_1-lg1-lib1\tSM:H_UJ-GPD1001_NTN001-M1502596_1'
#Some of this data was never imported into the GMS and will be imported from a LIMS spreadsheet (obtained by work order)
#Some of this data is coming from the GMS and a GMS query will be used to obtain it
&getIdInfo('-sample_data'=>$sample_data, '-lims_data_file'=>$lims_data_file, '-assets'=>\%assets);

#ALIGNMENT JOBS - all samples
#Create an output directory for each alignment job
&createAlignDirs('-base_dir'=>$base_dir, '-alignment_dir'=>$alignment_dir, '-sample_data'=>$sample_data);

#Create a YAML file for each alignemt job
&createAlignYmls('-base'=>$alignment_dir, '-sample_data'=>$sample_data, '-yml_name'=>$alignment_yml_name,
                 '-reference'=>$reference, '-dbsnp'=>$dbsnp, '-known_indels'=>$known_indels, '-mills'=>$mills);

#Create a Toil command for each alignment job
&createAlignToils('-sample_data'=>$sample_data, '-lsb_additional'=>"docker(mgibio/unaligned-to-bqsr)", '-assets'=>\%assets, 
                  '-align_cwl_exome'=>$align_cwl_exome, '-align_cwl_swift'=>$align_cwl_swift, 
                  '-align_bash_file'=>$align_bash_file, '-align_session_file'=>$align_session_file,
                  '-job_group'=>$job_group);


#SOMATIC JOBS - each tumor/normal sample pair for only the samples with exome data
#Create an output directory for each somatic job


#Create a YAML file for each somatic job 

#Create a Toil command for each somatic job

#Print out a TSV file that summarizes results labels and storage locations for results
&printLegend('-legend_file'=>$sample_legend_file, '-sample_data'=>$sample_data);


#Check for missing files and directories (i.e. all %assets)
foreach my $a (sort keys %assets){
  if ($assets{$a} eq 'f'){
    die "\n\nNeccessary file not found: $a\n\n" unless (-e $a);
  }
  if ($assets{$a} eq 'd'){
    die "\n\nNeccessary directory not found: $a" unless (-e $a && -d $a)
  }
}


print "\n\n";
exit;

sub loadData{
  my %args = @_;
  my $file = $args{'-data_file'};
  my $assets = $args{'-assets'};

  my %samples;
  my %unique_labels;
  open(INFILE, "$file") || die "\n\ncould not open infile: $file\n\n";
  my $o = 0;
  while(<INFILE>){
    next if $_ =~ /^Instrument/;
    chomp $_;
    my @l = split("\t", $_);
    my $id = $l[0];
    my $data_type = $l[1];
    my $i_name = $l[3];
    my $s_name = $l[4];
    my $label = $l[5];
    my $sc_name = $l[7];
    my $tissue = $l[9];
    my $e_type = $l[10];
    my $bam_path = $l[13];
    my $unique_label = $i_name . "_" . $label . "_" . $sc_name . "_" . $tissue . "_" . $e_type . "_" . $data_type;
    $unique_labels{$unique_label} = $s_name;
    my $s = $unique_label;

    if (defined($samples{$s})){
      my $data = $samples{$s}{data};
      $data->{$id}->{path} = $bam_path;
    }else{
      $o++;
      my %data;
      $data{$id}{path} = $bam_path;
      $samples{$s}{data} = \%data;
      $samples{$s}{sample_name} = $s_name;
      $samples{$s}{data_type} = $data_type;
      $samples{$s}{individual_name} = $i_name;
      $samples{$s}{sample_name} = $s_name;
      $samples{$s}{label} = $label;
      $samples{$s}{sample_common_name} = $sc_name;
      $samples{$s}{tissue_description} = $tissue;
      $samples{$s}{extraction_type} = $e_type;
      $samples{$s}{unique_label} = $unique_label;
      $samples{$s}{order} = $o;
    }

    #Make sure unique labels are unique
    die "\n\nUnique labels constructed are not unique\n\n" unless (keys %unique_labels == keys %samples); 

  }

  close(INFILE);

  return (\%samples);
}

sub determinePairs{
  my %args = @_;
  my $sample_data = $args{'-sample_data'};

  my %pairs;
  foreach my $s (sort {$sample_data->{$a}->{order} <=> $sample_data->{$b}->{order}} keys %{$sample_data}){
    my $i_name = $sample_data->{$s}->{individual_name};
    my $data_type = $sample_data->{$s}->{data_type};
    my $sample_common_name = $sample_data->{$s}->{sample_common_name};
    next unless ($data_type eq 'Exome');
    if ($sample_common_name eq 'normal'){
      $pairs{$i_name}{normal} = $s;
    }
  }

  foreach my $s (sort {$sample_data->{$a}->{order} <=> $sample_data->{$b}->{order}} keys %{$sample_data}){
    my $i_name = $sample_data->{$s}->{individual_name};
    my $data_type = $sample_data->{$s}->{data_type};
    my $sample_common_name = $sample_data->{$s}->{sample_common_name};
    next unless ($data_type eq 'Exome');
    next if ($sample_common_name eq 'normal');
    if (defined($pairs{$i_name}{tumors})){
      my $tumors = $pairs{$i_name}{tumors};
      $tumors->{$s} = 1;
    }else{
      my %tumors;
      $tumors{$s} = 1;
      $pairs{$i_name}{tumors} = \%tumors;
    }
  }

  return(\%pairs);
}

sub printLegend{
  my %args = @_;
  my $file = $args{'-legend_file'};
  my $sd = $args{'-sample_data'};
  
  open(OUTFILE, ">$file") || die "could not open legend file: $file";
  print OUTFILE "sample_name\tunique_label\tdata_type\tindividual_name\tsample_name\tlabel\tsample_common_name\ttissue_description\textraction_type\tinstrument_data_count\n";

  foreach my $s (sort {$sd->{$a}->{order} <=> $sd->{$b}->{order}} keys %{$sd}){
    my $id_count = keys %{$sd->{$s}->{data}};
    print OUTFILE "$sd->{$s}->{sample_name}\t$sd->{$s}{unique_label}\t$sd->{$s}{data_type}\t$sd->{$s}{individual_name}\t$sd->{$s}{sample_name}\t$sd->{$s}{label}\t$sd->{$s}{sample_common_name}\t$sd->{$s}{tissue_description}\t$sd->{$s}{extraction_type}\t$id_count\n";
  }
  close(OUTFILE);

  return;
}

sub createAlignDirs{
  my %args = @_;
  my $base_dir = $args{'-base_dir'};
  my $alignment_dir = $args{'-alignment_dir'};
  my $sample_data = $args{'-sample_data'};

  foreach my $s (sort keys %{$sample_data}){
    my $unique_label = $sample_data->{$s}->{unique_label};
    my $dir1 = $alignment_dir . $unique_label;
    $sample_data->{$s}->{dir} = $dir1;
    my $dir2 = $base_dir . "/workDir/" . $unique_label;
    $sample_data->{$s}->{work_dir} = $dir2;

    #If the run is marked as complete, do nothing
    my $complete_file = $dir1 . "/COMPLETE";
    next if (-e $complete_file);

    my $cmd1 = "mkdir $dir1";
    unless (-e $dir1 && -d $dir1){
      print "\n$cmd1";
      system($cmd1);
    }
    
    my $cmd2 = "mkdir $dir2";
    unless (-e $dir2 && -d $dir2){
      print "\n$cmd2";
      system($cmd2);
    }
  }
  return;
}

sub getIdInfo{
  my %args = @_;
  my $sample_data = $args{'-sample_data'};
  my $lims_data_file = $args{'-lims_data_file'};
  my $assets = $args{'-assets'};

  use Genome;

  #Info needed:
  # '@RG\tID:2895804954\tPL:ILLUMINA\tPU:H75J7BBXX-CATTTTAT-CATTTTAT.8\tLB:H_UJ-GPD1001_NTN001-M1502596_1-lg1-lib1\tSM:H_UJ-GPD1001_NTN001-M1502596_1'

  #Load data manually extracted from the LIMS system
  my %lims_data;
  open (LD, "$lims_data_file") || die "\n\nCould not open lims data file: $lims_data_file\n\n";
  while(<LD>){
    next if ($_ =~ /^Instrument/);
    chomp $_;
    my @ld = split("\t", $_);
    my $id = $ld[0];
    $lims_data{$id}{flowcell} = $ld[1];
    $lims_data{$id}{lane} = $ld[2];
    $lims_data{$id}{index_sequence} = $ld[3];
    $lims_data{$id}{library_name} = $ld[6];
  }
  close(LD);

  #Go through each piece of instrument data and query GMS for info, if not found check the LIMS file, if neither works throw error
  foreach my $s (sort keys %{$sample_data}){
    my $data = $sample_data->{$s}->{data};
    foreach my $id (sort keys %{$data}){
      my $flowcell; 
      my $index_sequence;
      my $lane;
      my $library_name;

      my $instdata = Genome::InstrumentData->get($id);
      if ($lims_data{$id}){
        $flowcell = $lims_data{$id}{flowcell};
        $index_sequence = $lims_data{$id}{index_sequence};
        $lane = $lims_data{$id}{lane};
        $library_name = $lims_data{$id}{library_name};
      }elsif($instdata){
        $flowcell =  $instdata->flow_cell_id();
        $index_sequence = $instdata->index_sequence();
        $lane = $instdata->lane();
        $library_name =  $instdata->library->name();
        $data->{$id}->{path} = $instdata->bam_path;
      }else{
        print "\nCould not find info for instrument data id: $id";
      }
      my $read_group_string = '@RG\tID:'. $id .'\tPL:ILLUMINA\tPU:'. $flowcell .'-'. $index_sequence .'.'. $lane .'\tLB:'. $library_name .'\tSM:'. $sample_data->{$s}->{sample_name};
      $data->{$id}->{read_group_string} = $read_group_string;
      $assets->{$data->{$id}->{path}} = "f";
    }
  }

  return;
}

sub createAlignYmls{
  my %args = @_;
  my $alignment_dir = $args{'-base'};
  my $sample_data = $args{'-sample_data'};
  my $yml_name = $args{'-yml_name'};
  my $reference = $args{'-reference'};
  my $dbsnp = $args{'-dbsnp'};
  my $known_indels = $args{'-known_indels'};
  my $mills = $args{'-mills'};

  foreach my $s (sort keys %{$sample_data}){
    my $unique_label = $sample_data->{$s}->{unique_label};
    my $dir = $alignment_dir . $unique_label . "/";
    
    my $yml_file = $dir . $yml_name;
    $sample_data->{$s}->{yml_path} = $yml_file;

    #If the run is marked as complete, do nothing
    my $complete_file = $dir . "/COMPLETE";
    next if (-e $complete_file);

    my $data = $sample_data->{$s}->{data};
    my @read_groups;
    my @bams;
    foreach my $id (sort keys %{$data}){
      my $bam_path = $data->{$id}->{path};
      my $read_group_string = $data->{$id}->{read_group_string};
      push (@read_groups, "\'$read_group_string\'");
      push (@bams, "    - {class: File, path: $bam_path}");
    }
    my $read_groups = join(",", @read_groups);
    my $bams = join("\n", @bams);


my $yml = <<EOF;
reference:
    class: File
    path: $reference
bams:
$bams
readgroups: [$read_groups]
dbsnp:
    class: File
    path: $dbsnp
known_indels:
    class: File
    path: $known_indels
mills:
    class: File
    path: $mills
EOF
 
    open (YML, ">$yml_file") || die "\n\nCould not open out file: $yml_file\n\n";
    print YML $yml;
    close(YML);
  }

  return;
}

sub createAlignToils{
  my %args = @_;
  my $sample_data = $args{'-sample_data'};
  my $lsb_additional = $args{'-lsb_additional'};
  my $assets = $args{'-assets'};
  my $align_cwl_exome = $args{'-align_cwl_exome'};
  my $align_cwl_swift = $args{'-align_cwl_swift'};
  my $align_bash_file = $args{'-align_bash_file'};
  my $align_session_file = $args{'-align_session_file'};
  my $job_group = $args{'-job_group'};

  open(CWLBASH, ">$align_bash_file") || die "\n\ncould not open align cwl bash output file: $align_bash_file\n\n";
  my $j = 0;
  foreach my $s (sort keys %{$sample_data}){
    my $align_yml = $sample_data->{$s}->{yml_path}; $assets->{$align_yml} = "f";
    my $dir = $sample_data->{$s}->{dir}; $assets->{$dir} = "d";
    my $complete_file = $dir . "/COMPLETE";
    my $work_dir = $sample_data->{$s}->{work_dir}; 
    $assets->{$work_dir} = "d" unless (-e $complete_file);
    my $align_cwl;
    if ($sample_data->{$s}->{data_type} eq 'Exome'){
      $align_cwl = $align_cwl_exome;
    }elsif($sample_data->{$s}->{data_type} eq 'Swift-Capture'){
      $align_cwl = $align_cwl_swift;
    }else{
      die "\n\nunrecognized data type: " . $sample_data->{$s}->{data_type}
    }
    $j++;
    my $job_name = "a_toil_" . $j;
    my $bsub_cmd = "bsub -J $job_name -g /mgtoil -o $dir/lsf.out -e $dir/lsf.err bash -c \"LSB_SUB_ADDITIONAL=\'$lsb_additional\' cwltoil --disableCaching --stats --logLevel=DEBUG --logFile=$dir/toil.log --outdir=$dir --workDir=$work_dir --jobStore=$base_dir/jobStore/$s --batchSystem lsf $align_cwl $align_yml\"";

    print CWLBASH $bsub_cmd . "\n";
  }
  close(CWLBASH);

  #Create notes for logging in with the neccessary session
  open(SESSION, ">$align_session_file") || die "\n\ncould not open session file: $align_session_file\n\n";
  print SESSION "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q docker-interactive -a \'docker(mgibio/unaligned-to-bqsr)\' -Is \$SHELL\n\n";
  print SESSION "export LSB_SUB_ADDITIONAL=\"docker(mgibio/unaligned-to-bqsr)\"\n";
  print SESSION "export LSB_DEFAULTQUEUE=\"research-hpc\"\n";
  print SESSION "export LSF_SERVERDIR=\"/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc\"\n";
  print SESSION "export LSF_LIBDIR=\"/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/lib\"\n";
  print SESSION "export LSF_BINDIR=\"/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin\"\n";
  print SESSION "export LSF_INCLUDEDIR=\"/opt/lsf9/9.1/include\"\n";
  print SESSION "export LSF_ENVDIR=\"/opt/lsf9/conf\"\n";
  print SESSION "export PATH=\"/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/etc:/opt/lsf9/9.1/linux2.6-glibc2.3-x86_64/bin:\$PATH\"\n";
  close(SESSION);

  return;
}





