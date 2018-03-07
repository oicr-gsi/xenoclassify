use strict;
use warnings;

(open my $FILE, "${ARGV[0]}");
(open my $OUTPUT, ">", "${ARGV[1]}");
my $count=1;

while (my $line=<$FILE>){
	chomp $line;
	if ($count%4==1) {
		my @fields=split / /,$line;
		print $OUTPUT "$fields[0]\n"
	}
	$count++;
}