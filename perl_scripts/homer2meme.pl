#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

my $input       = '';
my $extract     = '';
my $pseudocount = 0.01;

GetOptions(
    'i=s' => \$input,
    'e=s' => \$extract,
    'a=f' => \$pseudocount,
    'h'   => sub { usage() },
) or usage();

usage() unless $input;
die "Error: -a must be > 0.\n" unless $pseudocount > 0;

my $fh;
if ($input eq '-') {
    $fh = \*STDIN;
} elsif ($input =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($input)
        or die "Cannot open $input: $GunzipError";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!";
}

my $header_printed = 0;
my $in_motif = 0;
my $motif_id = '';
my $description = '';
my @matrix;

while (<$fh>) {
    chomp;
    next unless length($_);

    # HOMER header line starts with '>'
    if (/^>(.*)/) {
        my $rest = $1;

        # Flush previous motif
        if ($in_motif && @matrix) {
            print_meme_header() unless $header_printed;
            $header_printed = 1;
            print_meme_motif($motif_id, $description, \@matrix);
        }
        @matrix = ();

        my @parts = split /\t/, $rest;
        my $mid  = defined $parts[0] ? $parts[0] : 'motif';
        my $desc = defined $parts[1] ? $parts[1] : $mid;

        # Apply -e filter
        if ($extract && $mid ne $extract && $desc ne $extract) {
            $in_motif = 0;
            next;
        }

        $motif_id    = $mid;
        $description = $desc;
        $in_motif    = 1;
        next;
    }

    next unless $in_motif;

    # Matrix data row
    my @tokens = split /\s+/;
    my $all_numeric = 1;
    for my $t (@tokens) {
        unless ($t =~ /^-?[\d.]+([eE][+-]?\d+)?$/) {
            $all_numeric = 0;
            last;
        }
    }
    next unless $all_numeric && @tokens;

    my @row = map { $_ + 0 } @tokens;
    if (scalar(@row) != 4) {
        warn "Warning: skipping malformed matrix row (expected 4 cols, got "
             . scalar(@row) . "): $_\n";
        next;
    }

    # Auto-detect log-odds vs probability (prob rows sum to ~1.0)
    my $sum = 0;
    $sum += $_ for @row;
    if ($sum < 0.98 || $sum > 1.02) {
        @row = logodds_to_prob(\@row, $pseudocount);
    }
    push @matrix, \@row;
}

# Flush last motif
if ($in_motif && @matrix) {
    print_meme_header() unless $header_printed;
    print_meme_motif($motif_id, $description, \@matrix);
}

# Close handle correctly depending on type
if ($input ne '-') {
    if (ref $fh && $fh->isa('IO::Uncompress::Gunzip')) {
        $fh->close() or warn "Error closing gz file: $GunzipError";
    } else {
        close $fh or warn "Error closing file: $!";
    }
}

# ---------------------------------------------------------------------------

sub logodds_to_prob {
    my ($row_ref, $pc) = @_;
    my $background = 0.25;
    my @raw = map { 2 ** $_ * $background } @$row_ref;
    my $total = $pc * scalar(@raw);
    $total += $_ for @raw;
    return map { ($_ + $pc) / $total } @raw;
}

sub print_meme_header {
    print "MEME version 4\n";
    print "\n";
    print "ALPHABET= ACGT\n";
    print "\n";
    print "strands: + -\n";
    print "\n";
    print "Background letter frequencies\n";
    print "A 0.25 C 0.25 G 0.25 T 0.25\n";
    print "\n";
}

sub print_meme_motif {
    my ($id, $desc, $matrix_ref) = @_;
    my $width = scalar @$matrix_ref;
    print "MOTIF $id $desc\n";
    print "\n";
    print "letter-probability matrix: alength= 4 w= $width nsites= 20 E= 0\n";
    foreach my $row (@$matrix_ref) {
        print "  " . join("  ", map { sprintf("%.6f", $_) } @$row) . "\n";
    }
    print "\n";
}

sub usage {
    print <<EOF;
Usage: $0 -i <input_file> [OPTIONS]

Convert HOMER motif format to MEME format.

Options:
    -i <file>   Input HOMER motif file (or '-' for stdin, supports .gz)
    -e <string> Extract only specified motif by id or description
    -a <float>  Pseudocount for log-odds to probability conversion (default: 0.01)
    -h          Show this help

Examples:
    $0 -i results/motifs.homer > raw/motifs.meme
    $0 -i results/motifs.homer.gz > raw/motifs.meme
    $0 -i results/motifs.homer -e "CTCF/Jaspar"
    cat motifs.homer | $0 -i -

EOF
    exit 0;
}
