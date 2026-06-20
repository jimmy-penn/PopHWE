#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Time::HiRes qw(gettimeofday);
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);

my $script = $0;

if (@ARGV < 2) {
    print "\nCompute genotype counts for each ancestry group.\n\n";
    print "Usage: $script <in_file> <output_file>\n";
    print "\n";
    exit;
}

my $inFile = $ARGV[0];
my $outFile = $ARGV[1];

my @t1 = gettimeofday;

unless (-e $inFile) {
    print "\nERROR: didn't find file $inFile!\n\n";
    exit;
}


my @saveLines = ();
GetPopGenoCnts($inFile);

open OUTFILE, ">$outFile" or die "\nERROR: couldn't open $outFile for writing!\n";
foreach my $line (@saveLines) {
    print OUTFILE "$line\n";
}
close OUTFILE;

print "Results saved to $outFile\n\n";

my @t2 = gettimeofday;
my $time = GetTimeDifference(\@t1, \@t2);
print "Time used $time\n";


sub GetPopGenoCnts
{
    my ($file) = @_;

	open FILE, $file or die "Couldn't open $file\n";

    my $numSnps = 0;
    my $head = <FILE>;
    my ($num, $rs, $chr, $pop, $ref, $alt, @cnts) = ();
    chomp $head;
    ($num, $rs, $chr, $pop, $ref, $alt, @cnts) = split /\t/, $head;

    print "SNP $rs pop $pop al $ref/$alt\n";
    my $numPops = @cnts/3;
    my $line = "rs$rs";

    for my $i (0..$numPops-1) {
        my $cntRefHom = $cnts[3*$i];
        my $cntHet = $cnts[3*$i+1];
        my $cntAltHom = $cnts[3*$i+2];
        my $pop = "$pop.$i";
        $line .= "\t$pop";
    }
    push @saveLines, $line;

    while(<FILE>) {
        chomp;
        ($num, $rs, $chr, $pop, $ref, $alt, @cnts) = split /\t/, $_;
        $line = "rs$num";
        for my $i (0..$numPops-1) {
            my $cntRefHom = $cnts[3*$i];
            my $cntHet = $cnts[3*$i+1];
            my $cntAltHom = $cnts[3*$i+2];
            my $cnts = "$cntRefHom/$cntHet/$cntAltHom";
            $line .= "\t$cnts";
        }
        push @saveLines, $line;

        $numSnps++;
#        last if ($numSnps > 1000);
    }
    close FILE;

    print "Total $numSnps SNPs\n";
}

sub GetTimeDifference
{
    my ($t1, $t2) = @_;

    my $msg = "";

    my $t1sec = $$t1[0];
    my $t1us  = $$t1[1];
    my $t2sec = $$t2[0];
    my $t2us  = $$t2[1];

    my $ds = $t2sec - $t1sec;
    my $dus = $t2us - $t1us;

    if ($dus < 0) {
        $dus += 1000000;
        $ds -= 1;
    }

    if ($ds < 0 || ($ds == 0 && $dus < 0)) {
        print "Error: T1 > T2\n";
        return $msg;
    }

    my $dds = $ds % 60;
    $dds += $dus/1000000;

    if ($ds == 0) {
        $msg .= "$dus micro seconds";
        return $msg;
    }

    my $dh = int($ds/3600);
    $ds = $ds % 3600;
    my $dm = int($ds/60);

    $msg .= "$dh hour" if ($dh > 0);
    $msg .= "s" if ($dh > 1);
    $msg .= " $dm minute" if ($dh > 0 || $dm > 0);
    $msg .= "s" if ($dm > 1);
    $msg .= " $dds seconds";

    return $msg;
}
