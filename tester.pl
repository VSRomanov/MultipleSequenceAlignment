
use lib "lib";
use ProgressiveAlignment;

## to extract from FASTA format
use Bio::SeqIO;

use warnings;
use strict;

## Input:
my $gap_open_penalty = -2;
my $gap_extend_penalty = -1;

## scoring function (match vs mismatch):
sub score_sub {
     return ($_[0] eq $_[1]) ? 1 : -1;
}


## getting sequences from FASTA file
my $file = shift; 		## gets the first parameter from the BASH command line (filename here)
$file="my_fasta_test2.fa"; 	## FASTA file
my $seqIO_object = Bio::SeqIO->new(-file => $file,      -format => "fasta");

my %sequences;
while (my $seq_object = $seqIO_object->next_seq){
	## in case of "bad" caracters in the id
	$seq_object->{"primary_seq"}{"display_id"}=~s/[^a-zA-Z\d]//g;
	
	## in case of the empty gene names
	if ($seq_object->id() eq ""){$seq_object->id("sequence")}
	
	## in case of the same gene names, add a number
	my $h=$seq_object->id();
	for (my $i=2; exists($sequences{$h});$i++){$h=$seq_object->id()."_".$i}
	$seq_object->id($h);
	
	$sequences{$seq_object->id()}=$seq_object->seq;
}


## Constructing the object
my $object = ProgressiveAlignment->new(score_sub => \&score_sub, gap_open_penalty => $gap_open_penalty, gap_extend_penalty =>$gap_extend_penalty);


## to create the 'Output' directory if not exists:
if ( !-d "Output") {
	mkdir "Output" or die "Error creating directory: 'Output'";
}


my $alignment = $object->align_mult_seq(sequences =>\%sequences, out_tree_path =>"Output/my_tree.newick", out_distance_matrix_path=>"Output/my_matrix.phylip");

print (alignment_to_string($alignment));


## sub for printing the aligned strings into the file/console
sub alignment_to_string {
	my @input_alignment=@{$_[0]};
	my $alignment_strings;
	foreach my $i (0..scalar(@{$input_alignment[0]})-1) {
		my $string="";
		foreach my $j (0..scalar(@input_alignment)-1) {
			my $letter = $input_alignment[$j][$i];
			$string.=$letter;
		}
		
		$alignment_strings.=$string."\n"; ## to go to the next line during printing
	}
	return($alignment_strings);
}
