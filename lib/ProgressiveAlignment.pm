package ProgressiveAlignment;

use warnings;
use strict;

## to create the guide tree
use Bio::TreeIO;
use List::Util qw[min max sum];		## for 'min()' and 'max()' functions
use File::Basename;
use Cwd;
use Bio::Root::IO;
use File::Path;
use Try::Tiny;
use Getopt::Long;
use Bio::TreeIO::newick;

use lib "lib";
require NeedlemanWunschModified;

my %default_constructor_parameters=(gap_open_penalty=>3,  gap_extend_penalty =>2, score_sub => \&default_score_sub);

sub default_score_sub{return ($_[0] eq $_[1]) ? 1 : -1}

my $bin_dir=File::Basename::dirname(Cwd::abs_path($0))."/";
my $neighbor_default=$bin_dir."neighbor";

## Constructor

sub new {
    my $class = shift;	
    my %params=@_;
    my %self;
    for my $param (keys(%default_constructor_parameters)){
	    if (exists($params{$param})){$self{$param}=$params{$param}}
	    else {$self{$param}=$default_constructor_parameters{$param}}
    }
    return bless \%self, $class;  
}

sub align_mult_seq{
	my $self = shift;
	my %params=@_;
	if ((!exists($params{"sequences"})) || (scalar(keys(%{$params{"sequences"}}))<2)){die("At least two sequences must be provided...\n")}

	my %profiles;
		
	## transforming each sequence from a text scalar into the array
	for my $seq_name (keys(%{$params{"sequences"}})){

		my @seq = split("", $params{"sequences"}{$seq_name});
		my @seq_array_structure;
		
		foreach my $l (0..scalar(@seq)-1) {
			my $element = $seq[$l];
			push (@seq_array_structure, [$element])
		}
		
		## in case of "bad" caracters
		$seq_name=~s/[^a-zA-Z\d]//g;
		if ($seq_name eq ""){die("Empty sequence names are not allowed...\n")}
		$profiles{$seq_name}=\@seq_array_structure;
	} 
	
	my %pairs;
	my @keys=sort(keys(%profiles));
	foreach my $i (0..scalar(@keys)-2) {
		foreach my $j ($i+1..scalar(@keys)-1) {
			$pairs{make_pair_name($keys[$i],$keys[$j])}={profile1_idx => $keys[$i] , profile2_idx => $keys[$j]};
		}
	}
	
	## alignments for each pair
	for my $pair_name (keys(%pairs)) {
		my @seq1_array_structure = @{$profiles{$pairs{$pair_name}{"profile1_idx"}}};
		my @seq2_array_structure = @{$profiles{$pairs{$pair_name}{"profile2_idx"}}};
	
		## adding 'Alignment', 'Score' and 'Distance' keys to each pair-hash  
		($pairs{$pair_name}{"Alignment"},$pairs{$pair_name}{"Score"})= $self->align_profiles(\@seq1_array_structure,\@seq2_array_structure);        
		$pairs{$pair_name}{"Distance"}=distance($pairs{$pair_name}{"Alignment"}, scalar(@seq1_array_structure), scalar(@seq2_array_structure));
	}
	
	## to print the distance score into the "my_matrix" txt-file
	if ($params{"out_distance_matrix_path"} ne ""){
		open(my $FH, '>', $params{"out_distance_matrix_path"}) or die "Could not open file ".$params{"out_distance_matrix_path"} ."\n";
		print $FH scalar(keys(%profiles)), "\n"; ## the first line: number of sequences
		my @seq_names=keys(%profiles);
		for (my $i=0;$i<@seq_names;$i++){
			my @line=(make_phylip_name($seq_names[$i])); ## the first element of the array is the output of the 'make_phylip_name' function
			for (my $j=0;$j<@seq_names;$j++){
				if ($i==$j){push(@line,0)}
				else {push(@line,$pairs{make_pair_name($seq_names[$i],$seq_names[$j])}{"Distance"})}
			}
			print($FH join(" ",@line)."\n");
		}
		close $FH;
	}
	
	
	## Tree reconstruction from the newick file
	my $tree=UPGMA(\%pairs,\%profiles);
	
	%profiles=(%profiles,map{$_ => $pairs{$_}{"Alignment"}}(keys(%pairs)));
	my %scores=map{$_ => $pairs{$_}{"Score"}}(keys(%pairs));
	
	my $final_alignment=$self->alignment_of_tree_node($tree->get_root_node,\%profiles,\%scores); 
	if ($params{"out_tree_path"} ne ""){
		Bio::TreeIO->new(-file => ">".$params{"out_tree_path"}, -format => "newick")->write_tree($tree);
	}
	return($final_alignment);
}


## Tree reconstruction (UPGMA hierarchical clustering method)
sub UPGMA{
	my @seq_names=keys(%{$_[1]});
	my %seq_names_to_index;
	my $i=0;
	my @nodes=map{$seq_names_to_index{$_}=$i++;my $ret=Bio::Tree::Node->new(-id => $_);$ret->{"n"}=1;$ret->{"branch_length_sum"}=0;$ret}(@seq_names);
	my @pairs;
	for my $key(keys(%{$_[0]})){
		my %h=%{$_[0]{$key}};
		$h{"node_1"}=$nodes[$seq_names_to_index{$h{"profile1_idx"}}];
		$h{"node_2"}=$nodes[$seq_names_to_index{$h{"profile2_idx"}}];
		push(@pairs,\%h);
	}
	my $last_node;
	while (scalar(@pairs)>0){
		
		## to find a minimum distance among pairs
		my $min_index=undef;
		for (my $i=0; $i<@pairs;$i++){
			if((!defined($min_index)) || ($pairs[$i]{"Distance"}<$pairs[$min_index]{"Distance"})){$min_index=$i}
		}
		
		$last_node=Bio::Tree::Node->new();
		$last_node->add_Descendent($pairs[$min_index]{"node_1"});
		$pairs[$min_index]{"node_1"}->branch_length($pairs[$min_index]{"Distance"}/2-$pairs[$min_index]{"node_1"}{"branch_length_sum"});
		$last_node->add_Descendent($pairs[$min_index]{"node_2"});
		$pairs[$min_index]{"node_2"}->branch_length($pairs[$min_index]{"Distance"}/2-$pairs[$min_index]{"node_2"}{"branch_length_sum"});
		$last_node->{"branch_length_sum"}=$pairs[$min_index]{"Distance"}/2;
		$last_node->{"n"}=$pairs[$min_index]{"node_1"}->{"n"}+$pairs[$min_index]{"node_2"}->{"n"};
		push(@nodes,$last_node);
		
		## to find all the rest pairs except the one with the minimum distance
		my @new_pairs;
		my %new_distance_nodes;
		my @new_distance_nodes;
		for (my $i=0; $i<@pairs;$i++){
			my $b1=0;
			my $b2=0;
			if (($pairs[$min_index]{"node_1"}!=$pairs[$i]{"node_1"}) && ($pairs[$min_index]{"node_2"}!=$pairs[$i]{"node_1"})){
				$b1=1;
				if (!exists($new_distance_nodes{$pairs[$i]{"node_1"}})) {
					$new_distance_nodes{$pairs[$i]{"node_1"}}=0;
					push(@new_distance_nodes,$pairs[$i]{"node_1"})
				}
			}
			if (($pairs[$min_index]{"node_1"}!=$pairs[$i]{"node_2"}) && ($pairs[$min_index]{"node_2"}!=$pairs[$i]{"node_2"})){	
				$b2=1;
				if (!exists($new_distance_nodes{$pairs[$i]{"node_2"}})) {
					$new_distance_nodes{$pairs[$i]{"node_2"}}=0;
					push(@new_distance_nodes,$pairs[$i]{"node_2"})
				}
				if ($b1 && $b2) {push(@new_pairs,$pairs[$i])};
			}
		}
		
		for my $new_distance_node(@new_distance_nodes){
			my $old_dist_index_1=undef;
			my $old_dist_index_2=undef;
			for (my $i=0; $i<@pairs;$i++){
				if ((($pairs[$min_index]{"node_1"}==$pairs[$i]{"node_1"}) && ($new_distance_node==$pairs[$i]{"node_2"})) 
				|| (($pairs[$min_index]{"node_1"}==$pairs[$i]{"node_2"}) && ($new_distance_node==$pairs[$i]{"node_1"}))) {
					$old_dist_index_1=$i;
				}
				elsif ((($pairs[$min_index]{"node_2"}==$pairs[$i]{"node_1"}) && ($new_distance_node==$pairs[$i]{"node_2"})) 
				|| (($pairs[$min_index]{"node_2"}==$pairs[$i]{"node_2"}) && ($new_distance_node==$pairs[$i]{"node_1"}))) {
					$old_dist_index_2=$i;
				}
				
			}
			my $new_Distance=($pairs[$min_index]{"node_1"}->{"n"}*$pairs[$old_dist_index_1]{"Distance"} + $pairs[$min_index]{"node_2"}->{"n"}*$pairs[$old_dist_index_2]{"Distance"} )/($pairs[$min_index]{"node_1"}->{"n"}+$pairs[$min_index]{"node_2"}->{"n"});
			push(@new_pairs,{node_1 => $last_node, node_2=> $new_distance_node, Distance=>$new_Distance});
		}
		@pairs=@new_pairs;
	}
	my $tree=Bio::Tree::Tree->new($last_node);
	return($tree)
}


sub align_profiles{
	my $self = shift;
	my @seq1_array_structure=@{$_[0]};
	my @seq2_array_structure=@{$_[1]};
	my @current_pair_alignment=();

	## scoring function:
	local *multiple_score_sub = sub {
		my $score=0;
		my ($operation,$seq1_pos,$seq2_pos,$extra_info)=@_;
			if ($operation eq "substitution"){
				for (my $i=0; $i<scalar(@{$seq1_array_structure[$seq1_pos]});$i++){
					for (my $j=0; $j<scalar(@{$seq2_array_structure[$seq2_pos]});$j++){
						my $char_seq1=$seq1_array_structure[$seq1_pos][$i];
						my $char_seq2=$seq2_array_structure[$seq2_pos][$j];
						if (($char_seq1 ne "-") && ($char_seq2 ne "-")){
							$score+=&{$self->{"score_sub"}}($char_seq1,$char_seq2)
						} elsif (($char_seq1 eq "-") && ($char_seq2 ne "-")) {
							$score+=(($seq1_pos>0) && (($extra_info eq "up_left") || ($seq1_array_structure[$seq1_pos-1][$i] eq "-"))) ? $self->{"gap_extend_penalty"} : $self->{"gap_open_penalty"};
						} elsif (($char_seq2 eq "-") && ($char_seq1 ne "-")) {
							$score+=(($seq2_pos>0) && (($extra_info eq "up_left") || ($seq2_array_structure[$seq2_pos-1][$j] eq "-"))) ? $self->{"gap_extend_penalty"} : $self->{"gap_open_penalty"};
						}
					}
				}
			} 
			elsif ($operation eq "gap_open") {
				my ($seq_array_struct_with_gap,$seq_array_struct_without_gap,$pos_in_seq_with_gap,$pos_in_seq_without_gap) = ($extra_info eq "shift_in_seq1") ? (\@seq2_array_structure,\@seq1_array_structure,$seq2_pos,$seq1_pos) : (\@seq1_array_structure,\@seq2_array_structure,$seq1_pos,$seq2_pos);
				
				## if the simbol in 'seq_without_gap' is "-", then this var is 0, and the score is not changed (+0):
				my $lines_in_seq_without_gap=sum(map{($_ eq "-")?0:1}(@{$seq_array_struct_without_gap->[$pos_in_seq_without_gap]}));
				for (my $i=0; $i<scalar(@{$seq_array_struct_with_gap->[0]});$i++){
					$score+= ((($pos_in_seq_with_gap>-1) && ($seq_array_struct_with_gap->[$pos_in_seq_with_gap][$i] eq "-")) ? $self->{"gap_extend_penalty"} : $self->{"gap_open_penalty"})*$lines_in_seq_without_gap; #(($pos_in_seq_without_gap==0) && ($pos_in_seq_with_gap==0)) || 
				}	
			} elsif ($operation eq "gap_extend") {
				my ($seq_array_struct_with_gap,$seq_array_struct_without_gap,$pos_in_seq_with_gap,$pos_in_seq_without_gap) = ($extra_info eq "shift_in_seq1") ? (\@seq2_array_structure,\@seq1_array_structure,$seq2_pos,$seq1_pos) : (\@seq1_array_structure,\@seq2_array_structure,$seq1_pos,$seq2_pos);
								
				## if the simbol in 'seq_without_gap' is "-", then this var is 0, and the score is not changed (+0):
				my $lines_in_seq_without_gap=sum(map{($_ eq "-")?0:1}(@{$seq_array_struct_without_gap->[$pos_in_seq_without_gap]}));
				$score+=scalar(@{$seq_array_struct_with_gap->[0]})*$lines_in_seq_without_gap*$self->{"gap_extend_penalty"};
			}
			
			return($score);
	};
		## subroutines for the output visualisation
	local *on_align=sub  {
		my @current_column=(@{$seq1_array_structure[$_[0]]},@{$seq2_array_structure[$_[1]]});
		unshift(@current_pair_alignment,\@current_column); ## 'unshift' since we go backwards
	};
	local *on_shift_seq1=sub {
		my $number_of_gaps=scalar(@{$seq2_array_structure[0]});
		my @current_column=(@{$seq1_array_structure[$_[0]]},("-")x$number_of_gaps);
		unshift(@current_pair_alignment,\@current_column); ## 'unshift' since we go backwards
	};
	local *on_shift_seq2=sub {
		my $number_of_gaps=scalar(@{$seq1_array_structure[0]});
		my @current_column=(("-")x$number_of_gaps,@{$seq2_array_structure[$_[0]]});
		unshift(@current_pair_alignment,\@current_column); ## 'unshift' since we go backwards
	};
	
	############### Alternative alignments ###################
	## Called when there's more than one way to construct the optimal alignment, with 1 argument which is a hashref enumerating the possibilities.
	## The hash may contain the following keys: 
	## "align". If this key exists, the optimal alstringsignment may align two sequence items. The key's value is an arrayref with the positions of the paired items in \@a and \@b, respectively.
	## "shift_in_seq1". If this key exists, the optimal alignment may align an item of the first sequence with a gap in the second sequence. The key's value is the position of the item in \@a.
	## "shift_in_seq2". If this key exists, the optimal alignment may align a gap in the first sequence with an item of the second sequence. The key's value is the position of the item in \@b.
	
	## Simple solution: returns only the one (first) possible alignment
	local *on_select_align=sub {
		my $option=(sort(keys (%{$_[0]}))) [0]; ##returns only the first possible alignment. 
												##sort() function is to have a constant result of that first possible alignment. 
		my_select(@_,$option);
		return($option)
	};
	
	local *my_select=sub {
		my ($hash_ref,$option)=@_;
		if ($option eq "align") {
			on_align($hash_ref->{"align"}->[0],$hash_ref->{"align"}->[1])
		} elsif ($option eq "shift_in_seq1") {
			on_shift_seq1($hash_ref->{"shift_in_seq1"})
		} elsif ($option eq "shift_in_seq2") {
			on_shift_seq2($hash_ref->{"shift_in_seq2"})
		}
	};
	
	## Constructing the object
	my $matcher = NeedlemanWunschModified->new(\&multiple_score_sub, $self->{"gap_open_penalty"});
	$matcher->gap_open_penalty($self->{"gap_open_penalty"} );
	$matcher->gap_extend_penalty($self->{"gap_extend_penalty"});
		
	## scalar 'score' is assigned to the return value of the function 'align' of the object 'matcher'
	my $score = $matcher->align(
	           \@seq1_array_structure,
	           \@seq2_array_structure,
	           {   align  => \&on_align,
	               shift_in_seq1 => \&on_shift_seq1,
	               shift_in_seq2 => \&on_shift_seq2,
				   select_align => \&on_select_align
	           });

	return(\@current_pair_alignment,$score);
	
}

sub make_pair_name{
	return(join("_",sort(@_)));
}

## subroutine for the calculation of distances
sub distance {
	my ($alignment_ref,$length_seq1,$length_seq2)=@_;
	my $number_of_matches=0;
	my $comparison_column=0;
	
	for (my $i=0; $i<scalar(@{$alignment_ref});$i++) {
		if (($alignment_ref->[$i][0] ne "-") && ($alignment_ref->[$i][1] ne "-")){
	  		if ($alignment_ref->[$i][0] eq $alignment_ref->[$i][1]) {
				$number_of_matches += 1;
			}
			$comparison_column += 1;
		}
	}
	my $distance=1-$number_of_matches/max($length_seq1,$length_seq2);
}


## to get the name of the sequence in the phylip format
sub make_phylip_name{
	my $name = $_[0];
	## to have exactly 10 characters in the sequence names
	if (length($name) > 10) {
		$name=substr($name, 0, 4) . '**' . substr($name, length($name)-4, length($name));
	} else {
		my $num_spaces = 10 - length($name);
		$name.=(" ")x$num_spaces;
	}
	return($name)
}


sub alignment_of_tree_node{#maybe we find a better name
	my ($self,$node,$profiles,$scores)=@_;
	my @childs=$node->each_Descendent();
	my @leaf_childs;
	my @inner_node_childs;
	for my $child(@childs){
		if ($child->is_Leaf()) {push(@leaf_childs,$child)}
		else {push(@inner_node_childs,$child)}
	}
	my $profile_child_1;
	my $profile_child_2;
	if (scalar(@inner_node_childs)==2) {
		$profile_child_1=$self->alignment_of_tree_node($inner_node_childs[0],$profiles);
		$profile_child_2=$self->alignment_of_tree_node($inner_node_childs[1],$profiles);
	} elsif (scalar(@inner_node_childs)==1){
		$profile_child_1=$self->alignment_of_tree_node($inner_node_childs[0],$profiles);	
		$profile_child_2=$profiles->{$leaf_childs[0]->id()};
	} else {
		my $pair_name=make_pair_name($leaf_childs[0]->id(),$leaf_childs[1]->id());  ##getter function 'id'
		$node->id($scores->{$pair_name}); ## setter function 'id'
		return($profiles->{$pair_name});
	}
	my ($alignment,$score)=$self->align_profiles($profile_child_1,$profile_child_2);
	$node->id($score); ## setting the profile scores to the nodes 
	return($alignment);
}

return(1);
