#!/usr/bin/perl -w

use strict;

open OUT, "> screenlog.txt" or die "Can't write screenlog: $!\n";

print OUT "Plate\tWell\tFlag\tComment\n";

my $comment = "this is a comment";

my @plates = (1,6,10);

my @well_rows = qw(A B C D E F G H I J K L M N O P);

my $max_well = 24;

foreach my $plate (@plates){
	my $well = 1;
	while($well <= $max_well){
		my $formatted_well = $well;
		if($well < 10){
			$formatted_well = '0' . $well;
		}
		foreach my $well_row (@well_rows){
			print OUT "$plate\t$well_row$formatted_well\tNA\t$comment\n";			
		}
		$well ++;
	}
}

close OUT;


