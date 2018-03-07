use strict;
use warnings;

# check number of arguments passed is correct 
# my $num_args = $#ARGV + 1; 
# if ($num_args != 2) { 
#     print "\nUsage: python extractMappingQuality.pl bam_file read_names\n"; 
#     exit; 
# } 

# open bam files and read name files
my $mouseBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/mm10/sorted_6816873.bam";
my $humanBam="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/bwa/hg19/sorted_6816873.bam";
my $readHPath="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/sorted_host.fastq";
my $readGPath="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/sorted_graft.fastq";
my $readBPath="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/sorted_both.fastq";
my $readAPath="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/sorted_ambiguous.fastq";
my $readNPath="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/sorted_neither.fastq";
(open my $MOUSE,"samtools view -F 256  $mouseBam | ") || die "unable to open mouse bam: $mouseBam"; 
(open my $HUMAN,"samtools view -F 256 $humanBam |") || die "unable to open human bam: $humanBam"; 
(open my $readH,"$readHPath");
(open my $readG,"$readGPath");
(open my $readB,"$readBPath");
(open my $readA,"$readAPath");
(open my $readN,"$readNPath");

# files to write to 
my $outputH="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/aligned/host.sam";
my $outputG="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/aligned/graft.sam";
my $outputB="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/aligned/both.sam";
my $outputA="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/aligned/ambiguous.sam";
my $outputN="/.mounts/labs/gsiprojects/gsi/Xenograft-Classifier/data/seqware_prov_report/projects/AMLXP/xenome/6816873/aligned/neither.sam";
(open my $host,">","$outputH") || die "unable to open $outputH";
(open my $graft,">","$outputG") || die "unable to open $outputG";
(open my $both,">","$outputB") || die "unable to open $outputB"; 
(open my $ambiguous,">","$outputA") || die "unable to open $outputA"; 
(open my $neither,">","$outputN") || die "unable to open $outputN"; 

#check read name against bam file
my $hLine=<$readH>; chomp $hLine;
my $gLine=<$readG>; chomp $gLine;
my $aLine=<$readA>; chomp $aLine;
my $bLine=<$readB>; chomp $bLine;
my $nLine=<$readN>; chomp $nLine;

while(my $mouseLine=<$MOUSE>){

    my $humanLine=<$HUMAN>;
    chomp $mouseLine;
    chomp $humanLine;

    my @mfields=split /\t/,$mouseLine;
    my @hfields=split /\t/,$humanLine;

    if ($mfields[0] eq $hLine)
    {
        print "host\n";
        $hLine=printToFile($mouseLine,$humanLine,$host,$readH,$MOUSE,$HUMAN);
    }
    elsif ($mfields[0] eq $gLine)
    {
        print "graft\n";
        $gLine=printToFile($mouseLine,$humanLine,$graft,$readG,$MOUSE,$HUMAN);
    }
    elsif ($mfields[0] eq $aLine)
    {
        print "amb\n";
        $aLine=printToFile($mouseLine,$humanLine,$ambiguous,$readA,$MOUSE,$HUMAN);
    }
    elsif ($mfields[0] eq $bLine)
    {
        print "both\n";
        $bLine=printToFile($mouseLine,$humanLine,$both,$readB,$MOUSE,$HUMAN);
    }
    elsif ($mfields[0] eq $nLine)
    {
        print "neither\n";
        $nLine=printToFile($mouseLine,$humanLine,$neither,$readN,$MOUSE,$HUMAN);
    }
    else 
    {
        print "ERROR\nhost:$hLine\ngraft:$gLine\nboth:$bLine\nambiguous:$aLine\nneither:nLine\nbam:$mfields[0]\n";
    }
}

# to pass by reference, prepend arguments with backslash "\"
sub printToFile {
    my ($mouseLine1, $humanLine1, $outputFile, $readFile,$MOUSE, $HUMAN) = @_;
    my $mouseLine2;
    my $humanLine2;
    $mouseLine2=<$MOUSE>;
    $humanLine2=<$HUMAN>;
    chomp $mouseLine2;
    chomp $humanLine2;
    print $outputFile "$mouseLine1\n$mouseLine2\n";
    print $outputFile "$humanLine1\n$humanLine2\n\n";
    my $line=<$readFile>;
    if (defined($line)) {
        chomp $line;
        return $line;
    }
    else {
        return "End of File";
    }
}

# close files
close $MOUSE;
close $HUMAN;
close $readH;
close $readG;
close $readB;
close $readN;
close $readA;
close $host;
close $graft;
close $neither;
close $both;
close $ambiguous;
