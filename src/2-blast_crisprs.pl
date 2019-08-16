#!/usr/bin/perl
use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# GLOBALS - organism of interest and BLAST query file (required).
my ( $SPECIES, $F1, $F2, $OUT, $PAM, $regex );
my $usage = "Usage:
  $0 -s input_species -f1 wt_fasta -f2 snp_fasta -o outputfilename -pam pamSeq

Quick Help:
  -s    'hs' or 'dm'
  -f1    FASTA file of wild-type designs
  -f2    FASTA file of variant (SNP) designs
  -o     outputfilename
  -pam   PAM sequence";

# Check flags.
GetOptions (
    's=s'   => \$SPECIES,
    'f1=s'  => \$F1,
    'f2=s'  => \$F2,
    'o=s'   => \$OUT,
    'pam=s' => \$PAM,
    help    => sub { pod2usage($usage) }
) or pod2usage(2);
pod2usage($usage) and exit unless $SPECIES and $F1 and $F2 and $OUT and $PAM;

my $sequences = read_fasta( 'fasta_files/' . $SPECIES . '.fasta' );

if ($PAM eq '-NGG') {
    $regex = '[ATGC]{21}GG';
} elsif ($PAM eq '-NAG') {
    $regex = '[ATGC]{21}AG';
} else {
    say $PAM and exit;
}

my ( $filename ) = $OUT;
blast_crispr({
    db        => 'blast_dbs/' . $SPECIES, 
    query     => $F1,
    format    => 6,
    word_size => 10,
    output    => $filename . '-wt_blast.txt'
});

blast_crispr({
    db        => 'blast_dbs/' . $SPECIES,
    query     => $F2,
    format    => 6,
    word_size => 10,
    output    => $filename . '-snp_blast.txt'
});

# Opens given $sequence_file and returns list of chromsome -> sequence.
sub read_fasta {
    my ( $sequence_file ) = @_;
    my %sequences;

    open( FASTA, '<', $sequence_file ) or die $!;
    local $/ = '>';
    my $useless = <FASTA>;
    while (<FASTA>) {
        chomp;
        if ( /(.*?)\n(.*)/s ) {
            my @header = split( ' ', $1 );
            my $chrom = $header[0];
            my $sequence = $2;
            $sequence =~ s/\s//g;
            $sequences{$chrom} = uc $sequence;
        }
    }
    close FASTA;
    return \%sequences;
}

# BLAST CRISPR and print results to report file.
sub blast_crispr {
    my ( $args ) = @_;

    # Sept 2016 BLAST v2.2 command:
    #     ./blast-2.2.26/bin/blastall -p blastn -d $db -i $query -m $format -W $wordSize 
    #    -a 4 -e 100 -G 5 -E 2 -F T -n T
    my $cmd = 'blastn -db ' . $args->{db}
        . ' -query ' . $args->{query}
        . ' -outfmt ' . $args->{format}
        . ' -word_size ' . $args->{word_size}
        . ' -num_threads 4 -evalue 100 -gapopen 5 -gapextend 2 -dust yes -task megablast 2>'
        . $filename . '-warnings.err';
    
    open( ALIGN, '>', $args->{output} ) or die $!;

    open( BLAST, "$cmd |" );
    while (<BLAST>) {
        chomp;
        next if /^#/;
        
        my (
            $query_seq, $chrom, $pident, $len, $mismatches, 
            $gaps, $q_start, $q_end, $s_start, $s_end
        ) = split;
        my ( $eight, $unique, $pam, $strand ) = ( 0, 0, 0, '+' );
        $chrom =~ s/lcl\|//g; # remove 'lcl|' from id

        my $padded_start = sprintf( '%08d', $s_start - $q_start + 1 );
        my $padded_end   = sprintf( '%08d', $s_end + 23 - $q_end );
        if ( $s_start > $s_end ) {
            $strand = '-';
            $padded_start = sprintf( '%08d', $s_start + $q_start - 1 );
            $padded_end   = sprintf( '%08d', $s_end - ( 23 - $q_end ) );
        }
        my $padded_qs = sprintf( '%02d', $q_start );
        my $padded_qe = sprintf( '%02d', $q_end );
        
        my ( $subject_seq, $prev_base ) = 
            get_subject_seq( $padded_end, $chrom, $padded_start );
        
        if ( defined $subject_seq and $subject_seq =~ /$regex/ ) {
            my $alignment = '';
            my @subject_pos = split( '', $subject_seq );
            my @query_pos   = split( '', $query_seq );
            for ( my $i = 0; $i < scalar @subject_pos; $i++ ) {
                $alignment .= $subject_pos[$i] eq $query_pos[$i] ? '|' : 'X';
                
                if ( $subject_pos[$i] ne $query_pos[$i] ) {
                    $eight++ if $i < 8;
                    $unique++ if ( $i >= 8 and $i < 20 );
                    $pam++ if $i >= 20;
                }
            }

            my $ot_type = 'Off-Target';
            if ( $eight + $unique <= 5 ) {
                if ( $eight + $unique + $pam == 0 ) {
                    $ot_type = 'On-Target';
                }
                say ALIGN join( "\t", 
                    $chrom, $padded_start, $padded_end, $padded_qs, $padded_qe, 
                    $query_seq, $strand, $prev_base, $subject_seq, $alignment,
                    $eight, $unique, $pam, $ot_type 
                );
            }
        }
    }
    close BLAST;
    close ALIGN;
}

# Returns alignment sequence ($subject) and base 1-nt upsteam ($prev_base) of DB.
sub get_subject_seq {
    my ( $end, $id, $start ) = @_;
    
    my ( $subject, $prev_base );
    if ( exists $sequences->{$id} ) {
        if ( $start < $end ) {
            $subject = substr( $sequences->{$id}, $start - 2, 24 );
        }
        else {
            $subject = substr( $sequences->{$id}, $end - 1, 24 );
            $subject = reverse $subject;
            $subject =~ tr/ACGT/TGCA/;
        }
        $prev_base = substr( $subject, 0, 1 );
        $subject   = substr( $subject, 1, 23 ); 
    }
    return ( $subject, $prev_base );
}
