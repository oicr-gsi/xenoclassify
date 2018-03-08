use strict; 
use warnings; 
 
(open my $file, $ARGV[0]) || die "Could not open $ARGV[0]"; 
(open my $outFile, ">", "extractComb-60-60-60-60_graft.txt"); 
my @mapCombination = (60,60,60,60); 
my @extractSet; 
my $lineNum = 1; 
my $count=1; 
my @set; 
my $empty;
while ($count < 51) 
{ 
    
    for (my $i=0; $i<4; $i++) 
    { 
            $set[$i] = <$file>; 
            chomp $set[$i]; 
            my @fields=split(' ', $set[$i]); 
            $extractSet[$i] = $fields[4]; 
    }        

    if ($extractSet[0] == $mapCombination[0] && $extractSet[1] == $mapCombination[1] && $extractSet[2] == $mapCombination[2] && $extractSet[3] == $mapCombination[3])  
    { 
            for (my $i=0; $i<4; $i++) 
            { 
                    print $outFile "$set[$i]\n"; 
            } 
            print $outFile "\n"; 
            $count++; 
    } 
    $empty=<$file>; 
} 