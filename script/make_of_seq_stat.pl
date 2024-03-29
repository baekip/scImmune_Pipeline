#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my ($in_config, $input_path);
#GetOptions(
#    'config=s' => \$in_config,
#    'input=s' => \$input_path
#);
#
#if (!defined $in_config or !-f $in_config){
#    printUsage();
#    die "ERROR! check your config file with -c option\n";
#}

GetOptions(
    'config=s' => \$in_config,
);

if (!defined $in_config){
    printUsage();
    die "ERROR! check your config file with -c option\n";
}

sub printUsage {
    print "perl $0 -c <config_file> \n";
}

print "Status: Making Sequencing Stat File \n";

#my $config_path = dirname (abs_path $in_config);
#$in_config="$config_path/$in_config";
my %info;
read_general_config ($in_config, \%info);

#############################################################
#0.preparation
#############################################################

my $project_path=$info{project_path};
my $result_path=$info{result_path};
my $delivery_id=$info{delivery_tbi_id};
my @delivery_list=split /\,/, $delivery_id;

#############################################################
#1.read fastq stats xls file 
#############################################################
my $rawdata_path = "$project_path/rawdata";
my $sequence_stat = "$rawdata_path/Sequencing_Statistics_Result.xls";
my $fastq_path = "$result_path/02_fastqc_orig";
#my $fastq_path="$project_path/result//02_fastqc_orig";


open my $fh_stat, '>', $sequence_stat or die;
print $fh_stat "TBI_ID\tDelivery_ID\tSub_Seq\tTotalReads\tTotalBases\tTotalBases(Gb)\tGC_Count\tGC_Rate\tN_ZeroReads\tN_ZeroReadsRate\tN5_LessReads\tN5_LessReadsRate\tN_Count\tN_Rate\tQ30_MoreBases\tQ30_MoreBasesRate\tQ20_MoreBases\tQ20_MoreBasesRate\n";

foreach my $id (@delivery_list) {
    my ($delivery_id,$tbi_id,$type_id) = split /\:/, $id;
    my $sample_path="$fastq_path/$tbi_id/";
    checkDir($sample_path);
    my @fastq_list = ('I1', 'R1', 'R2');
    
    foreach my $sub_id (@fastq_list){
        my $rst_xls = "$sample_path/$tbi_id\_$sub_id.fastq.gz.rst.xls";
        checkFile($rst_xls);
       
    #####Read column value#####
        open my $fh_rst, '<:encoding(UTF-8)', $rst_xls or die;
        my $header_line =  <$fh_rst>;
        chomp ($header_line);
        my @headers = split ('\t', $header_line);

        my %combined;
        while (my $row = <$fh_rst>){
            chomp $row;
            my @row_list = split /\t/, $row;
            for my $header (@headers){
                push (@{$combined{$header}}, shift(@row_list));
            }
        }
        
        my $total_reads = join ('', @{$combined{'#totalReadCnt'}});
        my $total_length = join ('', @{$combined{'totalLength'}});
        my $total_GC = join ('', @{$combined{'totalGCCnt'}}); 
        my $Nzero_reads = join ('', @{$combined{'NZeroIncludeReadCnt'}});
        my $N5_reads = join ('', @{$combined{'N5IncludeReadCnt'}}); 
        my $N_count = join ('', @{$combined{'totalNCnt'}});
        my $Q30_reads_R1 = join ('', @{$combined{'totalQ30ReadCnt'}}); 
        my $Q30_base_R1 = join ('', @{$combined{'totalQ30BaseCnt'}}); 
        my $Q20_reads_R1 = join ('', @{$combined{'totalQ20ReadCnt'}}); 
        my $Q20_base_R1 = join ('', @{$combined{'totalQ20BaseCnt'}}); 
    #####correct stats column#####
    
        my $sample_rawdata="$fastq_path/$tbi_id/$tbi_id\_$sub_id.fastq.gz";
        checkFile($sample_rawdata);
        my $fastq_header = `zcat $sample_rawdata | head -n 1 `;
        chomp ($fastq_header);
        my @header_list = split /\:/, $fastq_header;
        my $index=$header_list[-1];
        
        my $total_base_gb=changeGbp($total_length)." Gb";
        my $gc_rate=$total_GC/$total_length*100;$gc_rate=RoundXL($gc_rate,2)."%";
        my $Nzero_rate=$Nzero_reads/$total_reads*100;$Nzero_rate=RoundXL($Nzero_rate,2)."%";
        my $N5_rate=$N5_reads/$total_reads*100;$N5_rate=RoundXL($N5_rate,2)."%";
        my $N_rate=$N_count/$total_length*100;$N_rate=RoundXL($N_rate,2)."%";
        my $Q30_base=$Q30_base_R1;
        my $Q30_rate=$Q30_base/$total_length*100;$Q30_rate=RoundXL($Q30_rate,2)."%";
        my $Q20_base=$Q20_base_R1;
        my $Q20_rate=$Q20_base/$total_length*100;$Q20_rate=RoundXL($Q20_rate,2)."%";
        print $fh_stat "$tbi_id\t$delivery_id\t$sub_id\t$total_reads\t$total_length\t$total_base_gb\t$total_GC\t$gc_rate\t$Nzero_reads\t$Nzero_rate\t$N5_reads\t$N5_rate\t$N_count\t$N_rate\t$Q30_base\t$Q30_rate\t$Q20_base\t$Q20_rate\n";
    }
}
close ($fh_stat);
print "Status: Making Seuquencing Stat - End\n"; 
#############################################################
#sub
#############################################################
sub changeGbp{
    my $val = shift;
    $val = $val/1000000000;
    $val = &RoundXL($val,2);
    return $val;
}

sub RoundXL {
    sprintf("%.$_[1]f", $_[0]);
}

sub checkDir{
    my $dir=shift;
    if (!-d $dir){
        die "ERROR! not exist <$dir>\n";
    }
}

sub checkFile{
	my $file = shift;
	if (!-f $file){
		die "ERROR ! not found <$file>\n";
	}
}

sub read_general_config{
	my ($file, $hash_ref) = @_;
	open my $fh, '<:encoding(UTF-8)', $file or die;
	while (my $row = <$fh>) {
		chomp $row;
		if ($row =~ /^#/){ next; } # pass header line
		if (length($row) == 0){ next; }

		my ($key, $value) = split /\=/, $row;
		$key = trim($key);
		$value = trim($value);
		$hash_ref->{$key} = $value;
	}
	close($fh);	
}

sub trim {
	my @result = @_;

	foreach (@result) {
		s/^\s+//;
		s/\s+$//;
	}

	return wantarray ? @result : $result[0];
}

sub make_dir {
    my $dir_name = shift;
    if (!-d $dir_name){
        my $cmd_mkdir = "mkdir -p $dir_name";
        system ($cmd_mkdir);
    }
}

