#!/usr/bin/perl
use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use Pod::Usage;
use POSIX 'log10';
use Data::Dumper;

# GLOBALS.
my ( $B1, $B2, $filename, $OUTPUT ); 
my $usage = "Usage:
  $0 -b1 wt_blast -b2 snp_blast [-out output_file]

Quick Help:
  -b1	BLAST report of wild-type designs
  -b2	BLAST report of variant (SNP) designs
  -out	output file (default: final_summary.csv)";

# Check flags.
GetOptions (
	'b1=s'  => \$B1,
	'b2=s'  => \$B2,
	'out:s' => \$filename,
	help    => sub { pod2usage($usage) }
) or pod2usage(2);
pod2usage($usage) and exit unless $B1 and $B2;

$OUTPUT = 'tmp/results/' . $filename . '.csv';

open( REPORT, '<', 'tmp/' . $filename . '-snp_summary.csv' ) or die $!;
my $header = <REPORT>;
chomp $header;

open( FINAL, '>', $OUTPUT ) or die $!;
say FINAL join( ',', $header, 'wt_off_target' , 'wt_efficiency', 'variant_off_target', 'variant_efficiency' );

our $on_target_hits;
my $wt_otes  = calculate_ote( $B1 );
my $snp_otes = calculate_ote( $B2 );
my $pssm = load_pssm( 'src/p_values.txt' );
while (<REPORT>) {
	chomp;

	# 733     752     +       750     A->G    TCACCACGATCAAAAGGAACAAG TCACCACGATCAAAAGGGACAAG
	my @columns = split ',';
	my $wt      = $columns[7];
	my $snp     = $columns[8];

	if ( $on_target_hits->{$wt} // 0 <= 1 && $on_target_hits->{$snp} // 0 == 0 ) {
		print FINAL $_;
		print FINAL ',', ( $wt_otes->{$wt} || '0' );
		print FINAL ',', score_match( $wt, $pssm );
		print FINAL ',', ( $snp_otes->{$snp} || '0' );
		print FINAL ',', score_match( $snp, $pssm );
		print FINAL "\n";
	}
}
close REPORT;
close FINAL;

sub calculate_ote {
	my ( $input ) = @_;

	my $ote_numbers = 'tmp/' . $filename . '_ote_numbers.txt';

	# Create temp file of OT scores pre-processed.
	`cut -f6,11-12,14 $input > $ote_numbers`;

	my $otes;
	open( OT, '<', $ote_numbers ) or die $!;
	while (<OT>) {
		chomp;
		my ( $crispr, $ot1, $ot2, $type ) = split "\t";
		if ( $type eq 'On-Target' ) {
			$on_target_hits->{$crispr}++;
			next;
		}
		my $total = $ot1 + $ot2;
		if ( $total != 0 and $total < 3 ) {
			$total = 3;
		}
		$otes->{$crispr}{$total}++;
	}
	close OT;

	# Remove temp file.
	`rm $ote_numbers`;

	my $scores;
	foreach my $crispr ( keys %$otes ) {
		my $ote3 = $otes->{$crispr}{3} // 0;
		my $ote4 = $otes->{$crispr}{4} // 0;
		my $ote5 = $otes->{$crispr}{5} // 0;
		$scores->{$crispr} = $ote3 + ( $ote4 / 10 ) + ( $ote5 / 100 );
	}
	return $scores;
}

# Returns matrix values from given file.
sub load_pssm {
	my ( $file ) = @_;
	my $matrix;

	open( SCORES, '<', $file ) or die $!;
	while (<SCORES>) {
		chomp;
		my @columns = split "\t";
		my $base = $columns[0];

		for ( my $position = 1; $position < scalar @columns; $position++ ) {
			$matrix->{$base}{$position} = $columns[$position];
		}
	}
	close SCORES;

	return $matrix;
}

# Scores given CRISPR ($match) against the given PSSM ($matrix).
sub score_match {
	my ( $match, $matrix ) = @_;

	my $score = 1;
	my @bases = split( '', $match );
	my $position = 1;

	foreach my $base ( @bases ) {
		my $pos_score = $matrix->{$base}{$position++};
		$score *= $pos_score if defined $pos_score;
	}
	return sprintf( '%.2f', log10( $score ) * -1 );
}
