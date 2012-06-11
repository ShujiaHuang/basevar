# Author : Shujia Huang
# Date   : 2013-01-10

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my (@posFile, $filelist, $regionFile); 
GetOptions (

	"i=s" => \@posFile,
	"l=s" => \$filelist,
	"r=s" => \$regionFile,
);
LoadFileFromList($filelist, \@posFile) if (defined $filelist);
Usage() unless @posFile > 0 && $regionFile;

my %region;
my %index;

GetRegion($regionFile, \%region, \%index);

for my $fn (@posFile) {
    print STDERR ">>> loading $fn >>>\n";
    open (IN, ($fn =~ /\.gz$/)? "gzip -dc $fn |": $fn) or die "Cannot open file $fn: $!\n";
    while (<IN>) {
    # Input format: chr10   105186  .  T  C
        chomp;
        my @tmp   = split;
        my $refId = $tmp[0];
        next if !exists $index{$refId};

        my $flag      = 1;
        my $isOverlap = 0;
        for(my $i = $index{$refId}; $i < @{$region{$refId}}; ++$i) {

            next if ($tmp[1] > $region{$refId}->[$i][1]);
            last if ($tmp[1] < $region{$refId}->[$i][0]); 

            if ($flag) {
                $flag = 0;
                $index{$refId} = $i;
            }
            $isOverlap = 1;
            last;
        }

        if (not $isOverlap) {
            print join "\t", @tmp[0..7]; print "\n";
        }

    }
    close IN;
}
print STDERR "** All done **";

###################################################################
###################################################################
sub GetRegion {

	my ($file, $region, $index) = @_;
	
	open I, $file or die "Cannot open file $file : $!\n";
	while (<I>) {
	# Input format:
	# chr1    65298905      65432187
		chomp;
		my @tmp = split;
		push (@{ $$region{$tmp[0]} }, [@tmp[1,2]]);
	}
	close I;

	for my $refId (keys %$region) {
		@{$$region{$refId}} = sort {$a->[0] <=> $b->[0]} @{$$region{$refId}};
		$$index{$refId}     = 0;
	}
}

sub OverlapLength {
    my ($reg1,$reg2,$tar1,$tar2)= @_;

    my $length = 0;
    if ($tar1 >= $reg1 && $tar1 <= $reg2 && $tar2 >= $reg2) {

        $length = $reg2 - $tar1 + 1;
    } elsif ($tar1 >= $reg1 && $tar1 <= $reg2 && $tar2 <= $reg2) {

        $length=$tar2-$tar1 + 1;
    } elsif ($tar2 >= $reg1 && $tar2 <= $reg2 && $tar1 <= $reg1) {

        $length = $tar2 - $reg1 + 1;
    } else {

        $length = $reg2 - $reg1 + 1;
    }
    return $length;
}

################################################################
sub IsOvlap {
# 
    my ($start1, $end1, $start2, $end2) = @_;
    my $isOvlp = 0;

    $isOvlp = 1 if $end1 >= $start2 and $start1 <= $end2;

    return $isOvlp;
}

sub OvlapReg {
# Reture the overlap region.
    my ($start1, $end1, $start2, $end2) = @_;
    my ($start, $end) = (0, 0);

    if ($end1 >= $start2 and $end1 <= $end2 and $start1 <= $start2) {
        $start = $start2;
        $end   = $end1;
    } elsif ($end1 >= $start2 and $end1 <= $end2 and $start1 > $start2) {
        die "The end position is bigger than start position : $end1 > $start1\n" if ($start1 > $end1);
        $start = $start1;
        $end   = $end1;
    } elsif ($start1 <= $start2 and $end1 > $end2) {
        die "The end position is bigger than start position : $end2 > $start2\n" if ($start2 > $end2);
        $start = $start2;
        $end   = $end2;
    } elsif ($start1 <= $end2 and $end1 > $end2) {
        $start = $start1;
        $end   = $end2;
    } else {
        die "[ERROR] Not right.\n";
	}

    return ($start, $end);
}


sub LoadFileFromList {

    my ($filelist, $files) = @_;
    open I, $filelist or die "Cannot open this file $filelist: $!\n";
    while (<I>) {
        chomp;
        push @$files, $_;
    }
    close I;
}


sub Usage {

	die qq/
Version : 0.0.1 ( 2012-01-10 )
Author  : Shujia Huang
Last Modify : 2013-01-10

    Usage : perl $0 -i [position file] -r [region file] > output

    Options :

        -i  [str]   Position file. Should be sorted! (Find the positions which in the Region file "-r".)
        -l  [str]   Input (position) file list. Should be sorted!
        -r  [str]   Region file. Not required sorting.

/;

}





