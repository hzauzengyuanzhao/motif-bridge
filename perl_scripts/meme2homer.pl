#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

my $input       = '';
my $db          = 'NA';
my $motif_name  = '';
my $extract     = '';
my $bg          = 0.25;   # -b: background probability (uniform default)
my $t_offset    = 4;      # -t: threshold offset in log2 bits

GetOptions(
    'i=s' => \$input,
    'j=s' => \$db,
    'k=s' => \$motif_name,
    'e=s' => \$extract,
    'b=f' => \$bg,
    't=f' => \$t_offset,
    'h'   => sub { usage() },
) or usage();

usage() unless $input;
die "Error: -b must be in (0, 1].\n" unless $bg > 0 && $bg <= 1;

my $fh;
if ($input eq '-') {
    $fh = \*STDIN;
} elsif ($input =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($input)
        or die "Cannot open $input: $GunzipError";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!";
}

my $in_motif  = 0;
my $in_matrix = 0;   # true only after "letter-probability matrix:" line
my $motif_id  = '';
my $description = '';
my @matrix;

while (<$fh>) {
    chomp;

    # ---- New MOTIF block ----
    if (/^MOTIF\s+(\S+)(?:\s+(.*))?/) {
        if ($in_motif && @matrix) {
            my $score = calculate_score(\@matrix, $bg, $t_offset);
            print_motif($motif_id, $description, $score, \@matrix);
        }

        $in_motif  = 1;
        $in_matrix = 0;
        $motif_id  = $1;
        my $original_name = defined $2 ? $2 : $1;
        $original_name =~ s/\s+/ /g; # normalize internal whitespace
        $original_name =~ s/^\s+|\s+$//g;

        $description = $motif_name ? "$motif_name/$db" : "$original_name/$db";

        # -e filter: skip motifs that don't match
        if ($extract && $motif_id ne $extract && $original_name ne $extract) {
            $in_motif = 0;  # explicitly off so subsequent lines are skipped
            @matrix   = ();
            next;
        }

        @matrix = ();
        next;
    }

    next unless $in_motif;

    # ---- letter-probability matrix header ----
    if (/^letter-probability matrix:/) {
        $in_matrix = 1;
        next;
    }

    # ---- End-of-motif markers ----
    if (/^URL/) {
        next;
    }
    if (/^\/\//) {
        if (@matrix) {
            my $score = calculate_score(\@matrix, $bg, $t_offset);
            print_motif($motif_id, $description, $score, \@matrix);
        }
        $in_motif  = 0;
        $in_matrix = 0;
        @matrix    = ();
        next;
    }

    # ---- Matrix data rows (only after letter-probability matrix: header) ----
    if ($in_matrix && /^\s*[\d.]/) {
        s/^\s+//;
        my @row = split /\s+/;
        # Validate: MEME probability rows must have exactly 4 columns (A C G T)
        if (scalar(@row) == 4) {
            push @matrix, \@row;
        } else {
            warn "Warning: skipping malformed matrix row (expected 4 cols, got "
                 . scalar(@row) . "): $_\n";
        }
    }
}

# Flush last motif
if ($in_motif && @matrix) {
    my $score = calculate_score(\@matrix, $bg, $t_offset);
    print_motif($motif_id, $description, $score, \@matrix);
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
# Subroutines
# ---------------------------------------------------------------------------

sub calculate_score {
    my ($matrix_ref, $bg, $t_offset) = @_;
    my $score = 0;
    foreach my $row (@$matrix_ref) {
        my $max_p = 0;
        foreach my $p (@$row) {
            $max_p = $p if $p > $max_p;
        }
        $score += log2($max_p / $bg) if $max_p > 0;
    }
    $score -= $t_offset;
    $score = 0 if $score < 0;
    return sprintf('%.6f', $score);
}

sub log2 {
    my ($x) = @_;
    return log($x) / log(2);
}

sub print_motif {
    my ($id, $desc, $score, $matrix_ref) = @_;
    # HOMER format: >id \t description \t threshold \t log-p-value \t pseudo-counts \t num-sites
    # 6 tab-separated fields are required; log-p, pseudo, sites default to 0.
    print ">$id\t$desc\t$score\t0\t0\t0\n";
    foreach my $row (@$matrix_ref) {
        print join("\t", map { sprintf('%.6f', $_) } @$row) . "\n";
    }
}

sub usage {
    print <<EOF;
Usage: $0 -i <input_file> [OPTIONS]

Convert MEME format to HOMER motif format.

Options:
    -i <file>    Input MEME format file (or '-' for stdin, supports .gz)
    -j <string>  Database name (default: NA)
    -k <string>  Motif name to use (default: name from MEME file)
    -e <string>  Extract only specified motif by id or name
    -b <float>   Background probability (default: 0.25, uniform)
    -t <float>   Threshold offset in log2 bits (default: 4)
    -h           Show this help

Examples:
    $0 -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer
    $0 -i raw/motifs.meme.gz -e MA0021.1
    cat motifs.meme | $0 -i -
    $0 -i motifs.meme -b 0.25 -t 6

EOF
    exit 0;
}
