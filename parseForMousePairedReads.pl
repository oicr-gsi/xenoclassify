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
    print "\nUsage: parseForMousePairedReads.pl Project SWID\n"; 
    exit; 
} 
 
# assign arguments to script variables 
# example input: 6816873 
my $mouseBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/mm10/${ARGV[1]}.bam"; 
my $humanBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/hg19/${ARGV[1]}.bam"; 
my $outputFile="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/bams/${ARGV[0]}/hg19/notMouse/twoReads/${ARGV[1]}.sam"; # change to bam if possible 
 
(open my $MOUSE,"samtools view -F 256  $mouseBam |") || die "unable to open mouse bam: $mouseBam"; 
(open my $HUMAN,"samtools view -F 256 $humanBam |") || die "unable to open human bam: $humanBam"; 
 
(open my $HUMANHEADER,"samtools view -H $humanBam |") || die "unable to open human header"; 
my @header=<$HUMANHEADER>; 
close $HUMANHEADER; 
 
 
(open my $OUT,">","$outputFile") || die "unable to open $outputFile"; 
print $OUT join("",@header); 
 
 
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

        $humanrec1=<$HUMAN>;
        $humanrec2=<$HUMAN>;
        chomp $humanrec1;
        chomp $humanrec2;
        my @hfields1=split /\t/,$humanrec1;
        my @hfields2=split /\t/,$humanrec2;

        # check that query names are the same
        if (!((($mfields1[0] eq $mfields2[0]) && ($hfields1[0] eq $hfields2[0])) && $mfields1[0] eq $hfields1[0])){
            die "record mismatch :\n$mouserec1\n$humanrec1\n$mouserec2\n$humanrec2";
        }
        
        ###### distinguish between reads that have not aligned, partially aligned, and fully aligned
        if ((($mfields1[1] & 4)==4) && (($mfields1[1] & 8)==8)){
            $mmapped=0; # not aligned
        }
        elsif ((($mfields1[1] & 4)==4) || (($mfields1[1] & 8)==8)){
            $mmapped=1; # one read aligned
        }
        else { 
            $mmapped=2; # both reads aligned
        }

        if ((($hfields1[1] & 4)==4) && (($hfields1[1] & 8)==8)){
            $hmapped=0; # not aligned
        }
        elsif ((($hfields1[1] & 4)==4) || (($hfields1[1] & 8)==8)){
            $hmapped=1; # one read aligned
        }
        else { 
            $hmapped=2; # both reads aligned
        }

        ###### classify segement as graft, host, ambiguous, both, or neither
        # both (aligns to mouse and human)
        if ($hmapped==2 && $mmapped==2){
            $bothCount++;
        }
        # neither (does not align to mouse or human)
        elsif ($hmapped==0 && $mmapped==0){
            $neitherCount++;
        }
        # host (only aligns to mouse, partially mouse and not human)
        elsif ($hmapped==0){
            $hostCount++;
        }
        # graft (only aligns to human, partially human and not mouse)
        elsif ($mmapped==0){
            $graftCount++;
        }
        # ambiguous (partially aligns to mouse and human, aligns mouse and partial human, aligns human and partial mouse)
        elsif ($hmapped==1 || $mmapped==1) {
            $ambiguousCount++;
        }

        ###### print non-mouse only reads to sam file
        if (!($hmapped==0))
        {
            print $OUT "$humanrec1\n";
            print $OUT "$humanrec2\n";
        }

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
close $OUT;
