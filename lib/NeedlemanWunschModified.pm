package NeedlemanWunschModified;

use warnings;
use strict;

use List::Util qw(max); 	## for 'max()' function
use Carp;		## for 'croak' function. Act like die() or warn(), but with a more useful message

use feature "say";


my $from_diag = 1;
my $from_up = 2;
my $from_left = 4;
my $from_diag_idx = 0;
my $from_up_idx = 1;
my $from_left_idx = 2;


## variables needed for the right calculation of the penalties in multiple sequence alignment (to distiguesh if it is 'indel' or 'substitution')
my $gap_open="gap_open";
my $gap_extend="gap_extend";
my $substitution="substitution";

## Constructor (the object is a pair of sequences)
sub new {
    my $class = shift;
    my $score_sub = shift;

    my $self = { score_sub => $score_sub}; 
    if (@_) {
        $self->{gap_penalty} = $_[0]; ## gap_penalty from the user
    }

    return bless $self, $class;
}

## attributing 'gap_open_penalty' value to the objectnew
sub gap_open_penalty {
    my $self = shift;

    if (@_) {
        $self->{gap_open_penalty} = $_[0];
    }

    return $self->{gap_open_penalty};
}

## attributing 'gap_extend_penalty' value to the object
sub gap_extend_penalty {
    my $self = shift;

    if (@_) {
        $self->{gap_extend_penalty} = $_[0];
    }

    return $self->{gap_extend_penalty};
}



## to check if all callbacks exist
sub _canonicalize_callbacks {
    my $callback;
    if (@_) {
        $callback = $_[0];
    } else {
        $callback = { };
    }

	## checks that all keys for the 'select_align' callback exist. If not, starts '_curry_callback' sub.
    if (exists($callback->{select_align})) {
        my @cn = qw(align shift_in_seq1 shift_in_seq2); ## the same as line 'my @cn=("align", "shift_in_seq1", "shift_in_seq2");'
		foreach (@cn) {
		    if (!exists($callback->{$_})) {
		        $callback->{$_} = _curry_callback($callback->{select_align}, $_);
		    }
		}
    }

    return $callback;
}

my $gap_open_end_gaps=0;
my $cheap_end_gaps=1;
my $free_end_gaps=2;

my $end_gap_option=$gap_open_end_gaps;
## Adding missing keys to the 'select_align' callback.
sub _curry_callback {
    my ($univ_callback, $missing_callback_name) = @_;
    my $callback;
    if ($missing_callback_name eq 'align') {

		## the following sub is a wrapper for the select_align method to create the missing callback ('align' here)
        $callback = sub {
		    my $arg = { align => [ @_ ] };
		    my $rv = &$univ_callback($arg);
		    croak "select_align callback returned invalid selection $rv."
		        unless $rv eq 'align';
    	};
    } else {
    	
		## the following sub is a wrapper for the select_align method to create the missing callback ('shift_in_seq1' or 'shift_in_seq2' here)
        $callback = sub {
		    my $arg = { $missing_callback_name => $_[0] };
		    my $rv = &$univ_callback($arg);
		    croak "select_align callback returned invalid selection $rv."
		        unless $rv eq $missing_callback_name;
		};
    }

    return $callback;
}


sub align {
    my $self = shift;

    my $seq1_ref = shift;
    my $seq2_ref = shift;

    $self->{callbacks} = _canonicalize_callbacks(@_); ## the key "callbacks" with the return of that sub as a value is added 
    												  ## to the $self (the object $matcher) 
    												  
    												  

	if ($self->{gap_open_penalty} >= $self->{gap_extend_penalty}) {
		warn "gap_open_penalty is smaller than gap_extend_penalty. Did you do this by purpose?";
    }
	return $self->_align_affine($seq1_ref, $seq2_ref);
        
}



## Alignment algorithm with affine gap penalty (two types of gap penalty: gap opening penalty and gap extension penalty)
## To find the best alignment, three similarity matrices are constructed by the scoring formula: 
##		 X(i,j) = {	max (D(i,j), U(i,j), L(i,j)) if i>0 & j>0;
##					L(i,j) if i=0 & j>0;
##					U(i,j) if i>0 & j=0;
## 					0 else	},
## where 
## 		"diag." matrix	D(i,j) = max (D(i-1,j-1), U(i-1,j-1), L(i-1,j-1)) + s(x_i,x_j),
## 		"up" matrix 	U(i,j) = max (D(i-1,j)+d, U(i-1,j)+e, L(i-1,j) + d),
## 		"left" matrix 	L(i,j) = max (D(i,j-1)+d, U(i,j-1)+d, L(i,j-1) + e),
##		s(x_i,x_j) - match/mismatch score at the (x_i,x_j) position,		
##		d - gap opening penalty,
##		e - gap extension penalty.
sub _align_affine {
    my $self = shift; ## the object $matcher
    my $seq1_ref = shift;
    my $seq2_ref = shift;

	my @A = ([ [ 0 ] ], [ [ 0 ] ], [ [ 0 ] ]); ## three "arrow" matrices with 0 at the starting (i=0, j=0) position
	
	## three score matrices with 0 at the starting (i=0, j=0) position
    my @D = ([ [ 0 ] ], [ [ 0 ] ], [ [ 0 ] ]); # indexed by $from_*_idx   ##(diag., up, left)
    my $length_seq1 = scalar(@$seq1_ref); ## seq1 will be on the top of matrices
    my $length_seq2 = scalar(@$seq2_ref); ## seq2 will be on the left of matrices


    my $score_diag = sub {
        my ($i, $j) = @_;
		my @base = map { $_->[$i - 1]->[$j - 1]; } @D; ## 'map' applies @D (3 matrices) to the 'diagonal-1' position (cell of the matrix)
		
		## when there is a gap caracter in the current column of the profile, we have to treat it as 'gap_extend' in case of that was up/left operation before or, in case of diagonal operation before, it depends on the column before (dealt within 'multiple_score_sub')  
		my @substitution_costs=(&{$self->{score_sub}}($substitution, $j - 1, $i - 1, "diagonal"),&{$self->{score_sub}}($substitution, $j - 1, $i - 1, "up_left"));  ## to distinguish 'substitution_costs' depending on the previous move (diagonal or up/left)
		@base=($base[0]+$substitution_costs[0],$base[1]+$substitution_costs[1],$base[2]+$substitution_costs[1]);
		my @which_max=which_max(@base);
		return($base[$which_max[0]],[$i-1,$j-1,\@which_max]);  ## first element - score for D, second element - arrow for A
    };

    my $score_up = sub {
		my ($i, $j) = @_;
		my @base = map { $_->[$i - 1]->[$j]; } @D;
		$base[$from_diag_idx] += &{$self->{score_sub}}($gap_open,$j-1,$i-1,"shift_in_seq2"); ## for matrix 1, the first element in the array @base
		$base[$from_up_idx] += &{$self->{score_sub}}($gap_extend,$j-1,$i-1,"shift_in_seq2"); ## for matrix 2, the second element in the array @base
		$base[$from_left_idx] += &{$self->{score_sub}}($gap_open,$j-1,$i-1,"shift_in_seq2"); ## for matrix 3, the third element in the array @base
		my @which_max=which_max(@base);
		return($base[$which_max[0]],[$i-1,$j,\@which_max]);  ## first element - score for D, second element - arrow for A
	};

    my $score_left;
	$score_left = sub {
	    my ($i, $j) = @_;
	    my @base = map { $_->[$i]->[$j - 1]; } @D;
	    $base[$from_diag_idx] += &{$self->{score_sub}}($gap_open,$j-1,$i-1,"shift_in_seq1");
	    $base[$from_up_idx] += &{$self->{score_sub}}($gap_open,$j-1,$i-1,"shift_in_seq1");
	    $base[$from_left_idx] += &{$self->{score_sub}}($gap_extend,$j-1,$i-1,"shift_in_seq1");
		my @which_max=which_max(@base);
		return($base[$which_max[0]],[$i,$j-1,\@which_max]);  ## first element - score for D, second element - arrow for A
	};
 

    ## filling the first rows of all three matrices
    my $j = 1; 
	while ($j <= $length_seq1) {
	    $D[$from_diag_idx]->[0]->[$j] = -"inf";
		$D[$from_up_idx]->[0]->[$j] = -"inf";
		$D[$from_left_idx]->[0]->[$j] = $D[$from_left_idx]->[0]->[$j-1] + &{$self->{score_sub}}((($j==1)?$gap_open:$gap_extend),$j-1,-1,"shift_in_seq1");
		$A[$from_diag_idx]->[0]->[$j] = undef;
		$A[$from_up_idx]->[0]->[$j] = undef;
		$A[$from_left_idx]->[0]->[$j] = [0,$j-1,[$from_left_idx]];
	    ++$j;
	}
     

    ## filling the first columns of all three matrices
    my $i = 1;
    while ($i <= $length_seq2) {
		$D[$from_diag_idx]->[$i]->[0] = -"inf";
		$D[$from_up_idx]->[$i]->[0] = $D[$from_up_idx]->[$i-1]->[0] + &{$self->{score_sub}}((($i==1)?$gap_open:$gap_extend),-1,$i-1,"shift_in_seq2");
		$D[$from_left_idx]->[$i]->[0] = -"inf";
		$A[$from_diag_idx]->[$i]->[0] = undef;
		$A[$from_up_idx]->[$i]->[0] = [$i-1,0,[$from_up_idx]];
		$A[$from_left_idx]->[$i]->[0] = undef;
		++$i;
    }

    ## order must correspond to $from_* constants
    my @subproblems = ( $score_diag, $score_up, $score_left );

    ## filling the matrices @D & A
    $i = 1;
    while ($i <= $length_seq2) {
		$j = 1;
		while ($j <= $length_seq1) {
		    my $k = 0;
		    while ($k < 3) { ## scalar(@D), scalar(@subproblems)
				($D[$k]->[$i]->[$j],$A[$k]->[$i]->[$j]) = &{$subproblems[$k]}($i, $j);
				++$k;
		    }
				
		    ++$j;
		}
		++$i;
    }

    my @three_scores = map { $_->[$length_seq2]->[$length_seq1]; } @D; ## 'map' applies @D (3 matrices) to the last positions (cells of the matrices)
    my @which_max=which_max(@three_scores);
    my $result_score = $three_scores[$which_max[0]];

    ## backtracking
    $i = $length_seq2;
    $j = $length_seq1;
    my $current_matrices=\@which_max;
    while (($i > 0) || ($j > 0)) {
    	
    	my @alt=map{[$A[$_][$i][$j][0],$A[$_][$i][$j][1]]}(@{$current_matrices});
	
		if (!@alt) {
		    die "internal error";
		}

		my $pointer = [ $i, $j ];  ## current position
		my $move;
		if (@alt == 1) { ## if the length of the array is one (only the one best alignment option)
		    $move = $self->_simple_trace_back($pointer, $alt[0],  ## '0' in the $alt[0], the position where we want to go, because length of @alt is 1
						      $self->{callbacks});
		} else {
		    $move = $self->_trace_back($pointer, \@alt);
		}
	
		if ($move eq 'align') {($i,$j,$current_matrices)=@{$A[0][$i][$j]};}
		elsif ($move eq 'shift_in_seq1') {($i,$j,$current_matrices)=@{$A[2][$i][$j]};}
		elsif ($move eq 'shift_in_seq2') {($i,$j,$current_matrices)=@{$A[1][$i][$j]};}

	}

    return $result_score;
}

## sub for backtracking when we have aternative alignments 
sub _trace_back {
    my ($self, $pointer, $alternative_next_coordinates) = @_;
    my $arg = {}; ## hashref of all possible moves with 'align', 'shift_in_seq1' and 'shift_in_seq2' as keys
    foreach my $next (@$alternative_next_coordinates) {
        my $move_option = $self->_simple_trace_back($pointer, $next, { });  ## empty '{ }' because of the structure
        																	## of arguments in the _simple_trace_back sub

		if ($move_option eq 'align') {
		    $arg->{align} = [ $pointer->[1] - 1, $pointer->[0] - 1 ]; ## adding this key 'align' to $arg
		} elsif ($move_option eq 'shift_in_seq1') {
		    $arg->{shift_in_seq1} = $pointer->[1] - 1; ## adding this key 'shift_in_seq1' to $arg
		} elsif ($move_option eq 'shift_in_seq2') {
		    $arg->{shift_in_seq2} = $pointer->[0] - 1; ## adding this key 'shift_in_seq2' to $arg
		} else {
		    die "internal error";
		}
    }

    my $move;
    my $callback = $self->{callbacks};
    if (exists($callback->{select_align})) {
        $move = &{$callback->{select_align}}($arg);
		if (!exists($arg->{$move})) {
		    die "select_align callback returned invalid selection $move.";
		}
    } else {
        my @cn = qw(align shift_in_seq1 shift_in_seq2);
		foreach my $move_option (@cn) {
		    if (exists($arg->{$move_option})) {
		        $move = $move_option;
			last;
		    }
	}

	if (!$move) {
	    die "internal error";
	}

	if (exists($callback->{$move})) {
	    if ($move eq 'align') {
	        &{$callback->{align}}(@{$arg->{align}});
	    } else {
	        &{$callback->{$move}}($arg->{$move});
	    }
	}
    }

    return $move;
}

## sub for backtracking when we have only one alignment
sub _simple_trace_back {
    my ($self, $pointer, $next, $callback) = @_;

    if ($next->[0] == $pointer->[0] - 1) {
        if ($next->[1] == $pointer->[1] - 1) {

			## this 'if' is working only in case of the length of the array @alt is one (only the one alignment option (see above))
		    if (exists($callback->{align})) {
		        &{$callback->{align}}($next->[1], $next->[0]);  ## Note: in the function, the first parameter refers to seq1,
		        												## the second parameter refers to seq2
		    }
	
		    return 'align';
		} else {
		    if ($next->[1] != $pointer->[1]) {
		        die "internal error";
		    }

			## this 'if' is working only in case of the length of the array @alt is one (only the one alignment option (see above)) 	
		    if (exists($callback->{shift_in_seq2})) {
		        &{$callback->{shift_in_seq2}}($pointer->[0] - 1);
	
		    }
	
		    return 'shift_in_seq2';
		}
    } else {
        if ($next->[0] != $pointer->[0]) {
	    	die "internal error";
		}

		if ($next->[1] != $pointer->[1] - 1) {
		    die "internal error";
		}

		## this 'if' is working only in case of the length of the array @alt is one (only the one alignment option (see above)) 
		if (exists($callback->{shift_in_seq1})) {
		    &{$callback->{shift_in_seq1}}($pointer->[1] - 1);
		}
	
		return 'shift_in_seq1';
    }
}

sub which_max{
	my @ret=(0);
	for (my $i=1; $i<@_;$i++){
		if ($_[$i]>$_[$ret[0]]){@ret=($i)}
		elsif($_[$i]==$_[$ret[0]]){push(@ret,$i)}
	}
	return(@ret);
}

1;
