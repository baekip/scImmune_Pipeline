#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;
use Sys::Hostname;
use File::Basename;
#use File::Basename qw(dirname);
use Cwd qw(abs_path);
use lib dirname(abs_path $0) . '/../library';
use Utils qw(read_config trim make_dir checkFile cmd_system);

my ($script, $program, $input_path, $sample, $sh_path, $output_path, $threads, $option, $config);
GetOptions (
    'script|s=s' => \$script,
    'program|p=s' => \$program,
    'input_path|i=s' => \$input_path,
    'sample|S=s' => \$sample,
    'log_path|l=s' => \$sh_path,
    'output_path|o=s' => \$output_path,
    'threads|t=s' => \$threads,
    'options|r=s' => \$option,
    'config_file|c=s' => \$config,
);

$input_path="$input_path/$sample";
my $hostname=hostname;

my %info;
read_config ($config, \%info);
my $fastp = $info{fastp};
my $script_path = $info{dev_path};
my $faster = "$script_path/util/FasterFastqStatistics";

if ( -d $output_path ) {
    my $cmd = "rm -r $output_path";
    system($cmd);
}

my ($option_orig, $q_value, $p_value) = split /\,/, $option;
make_dir ($sh_path);
make_dir ($output_path);

#TN1803L0441--TGCTGAGT-1_S8_L008_R1_001.fastq.gz > TN1803L0441_S8_L008_R1_001.fastq.gz

my @fastq_I1_list = glob ("$input_path/*_I1*.{fastq,fq}.gz");
my @fastq_R1_list = glob ("$input_path/*_R1*.{fastq,fq}.gz");
my @fastq_R2_list = glob ("$input_path/*_R2*.{fastq,fq}.gz");

my $unpaired_path = "$output_path/../unpaired/$sample/";
make_dir ($unpaired_path);

if (scalar (@fastq_R1_list) == 0 || scalar (@fastq_R2_list) == 0 || scalar(@fastq_I1_list) == 0){
    die "ERROR: Not exist!! check your cutadapt stat file <$input_path>";
}elsif (scalar (@fastq_R1_list) != scalar (@fastq_R2_list) || scalar (@fastq_R2_list) != scalar (@fastq_I1_list) ){
    die "ERROR: Not inconsistent!! check your cutadapt stat file I1, R1, R2 <$input_path>";
}else {
    
    foreach my $id (@fastq_R2_list){
#        my ($prefix, $ten_idx, $hiseq_idx, $type_idx, $tail) = split /\_/, $id;
#        my ($sample_id, $index) = split /--/, $prefix;
        my($filename, $directories, $suffix) = fileparse($id);
        my ($prefix, $ten_idx, $hiseq_idx, $type_idx, $tail) = split /\_/, $filename;
        my $I1_orig = "$prefix\_$ten_idx\_$hiseq_idx\_I1\_$tail";
        my $R1_orig = "$prefix\_$ten_idx\_$hiseq_idx\_R1\_$tail";
        my $R2_orig = "$prefix\_$ten_idx\_$hiseq_idx\_R2\_$tail";
        print "$id\n"; 
        print $I1_orig."\n";
        checkFile ("$input_path/$I1_orig");
        checkFile ("$input_path/$R1_orig");
        checkFile ("$input_path/$R2_orig");
        
        my $I1 = "$prefix\_$ten_idx\_$hiseq_idx\_I1\_$tail";
        my $R1 = "$prefix\_$ten_idx\_$hiseq_idx\_R1\_$tail";
        my $R2 = "$prefix\_$ten_idx\_$hiseq_idx\_R2\_$tail";
        
        my $sh_file = "$sh_path/fastp.$prefix\_$ten_idx.sh";
        open my $fh_sh, '>', $sh_file or die;
        print $fh_sh "#/bin/bash\n";
        print $fh_sh "#\$ -N fastp.$prefix\_$ten_idx\n";
        print $fh_sh "#\$ -wd $sh_path\n";
        print $fh_sh "#\$ -pe smp $threads\n";
        print $fh_sh "date\n";
        print $fh_sh "$fastp -q $q_value \\
        -u $p_value \\
        -w 4 \\
        --disable_adapter_trimming \\
        --disable_length_filtering \\
        -i $input_path/$R1_orig \\
        -o $output_path/$R1 \\
        -I $input_path/$R2_orig \\
        -O $output_path/$R2 \\
        -h $unpaired_path/$sample\.R1_R2.fastp.html\n\n";
        
        print $fh_sh "$fastp \\
        -w 4 \\
        --disable_adapter_trimming \\
        --disable_length_filtering \\
        --disable_quality_filtering \\
        -i $input_path/$I1_orig \\
        -o $output_path/$I1 \\
        -I $output_path/$R2 \\
        -O $unpaired_path/$R2 \\
        -h $unpaired_path/$sample\.I1_R2.fastp.html\n\n";
        
#        print $fh_sh "$faster $output_path/$I1\n";
#        print $fh_sh "$faster $output_path/$R1\n";
#        print $fh_sh "$faster $output_path/$R2\n";
        
        print $fh_sh "date\n";
        close $fh_sh;
        cmd_system ($sh_path, $hostname, $sh_file);
    }
} 

