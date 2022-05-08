
# Generic multiple sequence alignment

This Perl module enables multiple sequence alignment of any base objects such as symbols, numbers or data structures.
It is a considerably modified extension of [Algorithm::NeedlemanWunsch pairwise alignment
 module](https://metacpan.org/pod/Algorithm::NeedlemanWunsch) created by Vaclav Barta.

Match and mismatch scores, as well as affine gap penalty, can be defined freely, although there are defaults in the module.

Provided by Dr. Vasily S. Romanov and Dr. Arne Sahm as part of their work at the [Leibniz Institute on Aging](https://www.leibniz-fli.de/).

## Background

There are several algorithms for the purpose of aligning more than two sequences. However, progressive alignment is the most widely used approach. The distinctive feature of this implementation is that it offers progressive alignment in a generic way, i.e. sequence alignment of any base objects such as letters, numbers or even complex data structures is possible. An overrideable match/mismatch function is provided for this purpose.

[Progressive alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment#Progressive_alignment_construction) is conducted by combining pairwise alignments
beginning with the most similar pair and progressing to the most distantly related.
This process involves two main steps: 1) A guide tree is inferred to represent the similarity relationships between the sequences, and 2) the MSA is built by aligning sub-alignments sequentiallyfollowing the guide tree. In this implementation, the initial guide tree is determined by the [UPGMA clustering method](https://en.wikipedia.org/wiki/UPGMA).

Initial pair alignments as well as second-step pairwise additions to the constructed MSA are performed by two-sequence global alignment technique, the [Needleman–Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm). Briefly, at first,a similarity matrix is constructed using match and mismatch scores, as well as gap penalty, for each individual pair of base objects (e.g. symbols). 

Noteworthy, this implementation supports use of affine gap costs, i.e. two different gap costs/scores: "gap opening penalty" >= "gap extension penalty".

## Usage
Typical usage looks like this:


```
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
$file="my_fasta.fa"; 	## FASTA file
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
```

## Demo/Examples of output


Settings:\
&nbsp;&nbsp;&nbsp;&nbsp;$match = 1;\
&nbsp;&nbsp;&nbsp;&nbsp;$mismatch = -1;\
&nbsp;&nbsp;&nbsp;&nbsp;$gap_open_penalty = -2;\
&nbsp;&nbsp;&nbsp;&nbsp;$gap_extend_penalty = -1

---

&nbsp;\
&nbsp;\
Sequences to align 1:\
&nbsp;&nbsp;&nbsp;&nbsp;ATGTTGCACaACGCAGCCCT\
&nbsp;&nbsp;&nbsp;&nbsp;TgTtGCaCAACTCAGCCAtA\
&nbsp;&nbsp;&nbsp;&nbsp;ATGTGCAACGCAgCCcTA
```
ATGTTGCACaACGCAGCCCT-
ATG-TGCA--ACGCAgCCcTA
-TgTtGCaCAACTCAGCCAtA
```

&nbsp;\
&nbsp;\
Sequences to align 2:\
&nbsp;&nbsp;&nbsp;&nbsp;gene:GGCACCATAACACTTGTCACGTACTGG...protein:GTITLVTYW\
&nbsp;&nbsp;&nbsp;&nbsp;gene:GGCACGATAAGTTCTTTAGAAACGCAATGG...protein:GTISSLETQW\
&nbsp;&nbsp;&nbsp;&nbsp;gene:GGCACCAAAACGTTAAGGCCTACATACTGG...protein:GTKTLRPTYW
```
gene:GGCACCATAA---CACT--TGTC-ACGTACTGG...protein:GTI-TL-VTYW
gene:GGCACCAAAA---CGTTAAGGCCTACATACTGG...protein:GTK-TLRPTYW
gene:GGCACGATAAGTTCTTT--AGAA-ACGCAATGG...protein:GTISSL-ETQW
```

&nbsp;\
&nbsp;\
Sequences to align 3:\
&nbsp;&nbsp;&nbsp;&nbsp;1234567890\
&nbsp;&nbsp;&nbsp;&nbsp;2345678901\
&nbsp;&nbsp;&nbsp;&nbsp;3456789012\
&nbsp;&nbsp;&nbsp;&nbsp;4567x90123
```
1234567890---
-2345678901--
--3456789012-
---4567x90123
```

&nbsp;\
&nbsp;\
Sequences to align 4:\
&nbsp;&nbsp;&nbsp;&nbsp;1.Let's_see_if_it_works...\
&nbsp;&nbsp;&nbsp;&nbsp;2.Do_you_see_if_it_works_fine?..\
&nbsp;&nbsp;&nbsp;&nbsp;3.It's_fine!!!\
&nbsp;&nbsp;&nbsp;&nbsp;4.You_see?_It_is_fine!..
```
1.-Let's_see_if_it_works-----...
2.Do_you_see_if_it_works_fine?..
4.---You_see--?_It_---is_fine!..
3.--------------It'----s_fine!!!
```



