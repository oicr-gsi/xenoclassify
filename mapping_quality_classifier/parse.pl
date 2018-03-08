use strict;
use warnings;

# need to pass arguments to script 
# I could make this script a sub routine and pass arguments to 
# sub listy {
#    my ($foo) = @_;
# }
# Or I can pass arguments from command line to script 

# check number of arguments passed is correct
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: parse.pl Project SWID\n";
    exit;
}

# assign arguments to script variables
# example input: 6816873
my $mouseBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/mm10/${ARGV[1]}.bam";
my $humanBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/hg19/${ARGV[1]}.bam";
my $outputFile="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/hg19/notMouse/${ARGV[1]}.sam"; # change to bam if possible

(open my $MOUSE,"samtools view -F 256 $mouseBam |") || die "unable to open mouse bam: $mouseBam";
(open my $HUMAN,"samtools view -F 256 $humanBam |") || die "unable to open human bam: $humanBam";

(open my $HUMANHEADER,"samtools view -H $humanBam |") || die "unable to open human header";
my @header=<$HUMANHEADER>;
close $HUMANHEADER;


(open my $OUT,">","$outputFile") || die "unable to open $outputFile";
print $OUT join("",@header);


#my $reccount;
my $totalCount;
my $mouseCount;
my $humanCount;
my $bothCount;
my $neitherCount;

while(my $mouserec=<$MOUSE>){
	chomp $mouserec;
	my @mfields=split /\t/,$mouserec;

	my $humanrec=<$HUMAN>;
	chomp $humanrec;
	my @hfields=split /\t/,$humanrec;

	if($mfields[0] ne $hfields[0]){
		die "record mismatch :\n$mouserec\n$humanrec\n";
	}
	
	###### Bitwise operator &
	# As long as one pair of bits in the left and right operand are "1", the bitwise expression will evaluate to a number.
	# If this expression is in an if statement, it will also evaluate to a number, which evaluates to true.
	# Therefore, the following expression will evaluate to true as long as 0X4 or 0X8 is true.
	# Any reads with either of those flags will be marked as unmapped or "0".
	# This is a problem because the majority of reads will be marked as unmapped in mouse. The result is a low mouse only 
	# percentage (which is exactly what we see in the results).
	my $mmapped=$mfields[1] & 12 ? 0 : 1;
	my $hmapped=$hfields[1] & 12 ? 0 : 1;

	#print "$mfields[0] $mmapped $hmapped";<STDIN>;
	if($mmapped && !$hmapped){
		$mouseCount++;
	}
	else{
		print $OUT "$humanrec\n";
		if($mmapped && $hmapped)
		{
			$bothCount++;
		}
		elsif(!$mmapped && $hmapped) {
			$humanCount++;
		}
		else {
			$neitherCount++;
		}
	}

	$totalCount++;
}

my $percentageM = $mouseCount / $totalCount * 100;
my $percentageH = $humanCount / $totalCount * 100;
my $percentageB = $bothCount / $totalCount * 100;
my $percentageN = $neitherCount / $totalCount * 100;

print "Percentage of Mouse Reads: $percentageM\n";
print "Percentage of Human Reads: $percentageH\n";
print "Percentage of Both Reads: $percentageB\n";
print "Percentage of Neither Reads: $percentageN\n";

close $MOUSE;
close $HUMAN;
close $OUT; 




