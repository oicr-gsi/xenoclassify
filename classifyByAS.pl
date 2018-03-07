use strict;
use warnings;
 
# check number of arguments passed is correct 
my $num_args = $#ARGV + 1; 
if ($num_args != 2) { 
    print "\nUsage: classifyByAS.pl Project SWID\n"; 
    exit; 
} 
 
# assign arguments to script variables 
# example input: 6816873 
my $mouseBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/mm10/${ARGV[1]}.bam"; 
my $humanBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/${ARGV[1]}.bam"; 
# my $outputFile="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/hg19/notMouse/twoReads/${ARGV[1]}.sam"; # change to bam if possible 

# output files
(open my $GRAFT, ">", "/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/classifyByAS_results/${ARGV[1]}/graft.txt");
(open my $HOST, ">", "/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/classifyByAS_results/${ARGV[1]}/host.txt");
(open my $BOTH, ">", "/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/classifyByAS_results/${ARGV[1]}/both.txt");
(open my $NEITHER, ">", "/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/classifyByAS_results/${ARGV[1]}/neither.txt");
(open my $AMBIGUOUS, ">", "/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/${ARGV[0]}/bwa/hg19/classifyByAS_results/${ARGV[1]}/ambiguous.txt");
 
# open bams
(open my $MOUSE,"samtools view -F 256  $mouseBam |") || die "unable to open mouse bam: $mouseBam";  
(open my $HUMAN,"samtools view -F 256 $humanBam |") || die "unable to open human bam: $humanBam"; 
 
# (open my $HUMANHEADER,"samtools view -H $humanBam |") || die "unable to open human header"; 
# my @header=<$HUMANHEADER>; 
# close $HUMANHEADER; 
# (open my $OUT,">","$outputFile") || die "unable to open $outputFile"; 
# print $OUT join("",@header); 
 
# counters
my $totalCount=0; 
my $hostCount=0; 
my $graftCount=0; 
my $bothCount=0; 
my $neitherCount=0;
my $ambiguousCount=0; 

# read variables
my $mouserec1;
my $mouserec2;
my $humanrec1;
my $humanrec2;

# mapping type
my $mmapped;
my $hmapped;

while($mouserec1=<$MOUSE>){

        $mouserec2=<$MOUSE>;
        chomp $mouserec1;
        chomp $mouserec2;
        my @mfields1=split /\t/,$mouserec1;
        my @mfields2=split /\t/,$mouserec2;

        my $humanrec1=<$HUMAN>;
        my $humanrec2=<$HUMAN>;
        chomp $humanrec1;
        chomp $humanrec2;

        #get alignment scores
        my @AS_m1 = ($mouserec1 =~ m/AS:i:([0-9]*)/);
        my @AS_m2 = ($mouserec2 =~ m/AS:i:([0-9]*)/);
        my @AS_h1 = ($humanrec1 =~ m/AS:i:([0-9]*)/);
        my @AS_h2 = ($humanrec1 =~ m/AS:i:([0-9]*)/);

        my @hfields1=split /\t/,$humanrec1;
        my @hfields2=split /\t/,$humanrec2;

        # check that query names are the same
        if (!((($mfields1[0] eq $mfields2[0]) && ($hfields1[0] eq $hfields2[0])) && $mfields1[0] eq $hfields1[0])){
            die "record mismatch :\n$mouserec1\n$humanrec1\n$mouserec2\n$humanrec2";
        }
        
        ###### distinguish between reads that have not aligned, partially aligned, and fully aligned
        if (($AS_m1[0] > 104 || $AS_m2[0] > 104) && $AS_h1[0] < 104 && $AS_h2[0] < 104)
        {
            print $HOST "$hfields1[0]\n";
            $hostCount++
        }
        elsif ($AS_m1[0] < 104 && $AS_m2[0] < 104 && ($AS_h1[0] > 104 || $AS_h2[0] > 104))
        {
            print $GRAFT "$hfields1[0]\n";
            $graftCount++
        }
        elsif ($AS_m1[0] > 104 && $AS_m2[0] > 104 && $AS_h1[0] > 104 && $AS_h2[0] > 104)
        {
            print $BOTH "$hfields1[0]\n";
            $bothCount++
        }
        elsif ($AS_m1[0] < 10 && $AS_m2[0] < 10 && $AS_h1[0] < 10 && $AS_h2[0] < 10)
        {
            print $NEITHER "$hfields1[0]\n";
            $neitherCount++
        }
        else
        {
            print $AMBIGUOUS "$hfields1[0]\n";
            $ambiguousCount++
        }
        ###### print non-mouse only reads to sam file
        $totalCount++;
}

##### calculate class percentages
my $percentageH = ($hostCount / $totalCount) * 100;
my $percentageG = ($graftCount / $totalCount) * 100;
my $percentageB = ($bothCount / $totalCount) * 100;
my $percentageN = ($neitherCount / $totalCount) * 100;
my $percentageA = ($ambiguousCount / $totalCount) * 100;

$percentageH = sprintf("%.2f", $percentageH);
$percentageG = sprintf("%.2f", $percentageG);
$percentageB = sprintf("%.2f", $percentageB);
$percentageN = sprintf("%.2f", $percentageN);
$percentageA = sprintf("%.2f", $percentageA);

##### display results
print "Percentage of Mouse Reads: $percentageH\n";
print "Percentage of Human Reads: $percentageG\n";
print "Percentage of Both Reads: $percentageB\n";
print "Percentage of Neither Reads: $percentageN\n";
print "Percentage of Ambiguous Reads: $percentageA\n";
    
close $MOUSE;
close $HUMAN;
close $GRAFT;
close $HOST;
close $BOTH;
close $AMBIGUOUS;
close $NEITHER;
# close $OUT;
