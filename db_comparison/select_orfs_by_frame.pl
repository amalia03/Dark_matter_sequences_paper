#!/usr/bin/perl -w

## select sequences from the 6 frame translations that:
## 1. are split into: best negative / positive frames
##    That is select the best positive and negative frame.
## 2. ORFs that either go off either end of the reads; that
##    are consistent with the read containing either some
##    5' or 3' UTR.
## 3. That have a minimum length. If we only look at terminal
##    ORFs then we do not need to care so much about the coverage.
##    But do ouput the coverage into a report file of some sort.
## 4. Split output into different files.

($nuc_file, $pep_file, $min_length, $pos_odir, $neg_odir, $group_size) = @ARGV;

die "$pos_odir does not exist or is not a directory\n" if( ! -d $pos_odir );
die "$neg_odir does not exist or is not a directory\n" if( ! -d $neg_odir );

%nuc = read_fa( $nuc_file );
%pep = read_pep( $pep_file );

print_sequences( $pos_odir."/pep", $pep{fwd}, $group_size, $min_length );
print_sequences( $neg_odir."/pep", $pep{rev}, $group_size, $min_length );

## takes a vector of identifier lines. Returns 0
## on failuare and the id otherwise
sub check_ids {
    my @ids = @_;
    for(my $i=0; $i < @ids; $i++){
        $ids[$i] =~ s/^>(\S+).*$/$1/;
        return 0 if($ids[$i] ne $ids[0]);
    }
    return($ids[0]);
}

sub get_starts {
    my $id = shift @_;
    my $frame = shift @_;
    my @pep = (@_);
    die("no nucleotide sequence for $id\n"), if not defined($nuc{$id});
    my $seq = $nuc{$id};
    $seq = rev_complement($seq) if($frame < 0);
    $seq = substr( $seq, abs($frame) - 1 );
    my $pos = abs($frame);
    my @starts = ();
    for my $i(0..$#pep){
        push @starts, [ (substr($seq, 0, 3), $pos) ];
        my $offset = (1 + length($pep[$i])) * 3;
        $pos += $offset;
        $seq = substr($seq, $offset) if $offset < length($seq);
    }
    return(@starts);
}

sub extract_pep {
    my ($id, $seq, $frame) = @_;
    my @tmp = split /\*/, $seq;
    my @start_codons = get_starts($id, $frame, @tmp);
    return( [($tmp[0], 0, @{$start_codons[0]})] ) if(@tmp == 1);
    ## if more than one then return the longest orf.
    my @i = sort {length($tmp[$b]) <=> length($tmp[$a])} 0..$#tmp;
    my $i = $i[0];

 return( [($tmp[$i], -1, @{$start_codons[0]})] ) if $i == 0;
    return( [($tmp[$i], 1, @{$start_codons[0]})] ) if $i == $#tmp;
    return( [($tmp[$i], 0, @{$start_codons[0]})] );
}

sub extract_orfs {
    my @lines = @_;
    ## check the ids
    my $id = check_ids( @lines[0,2,4,6,8,10] );
    ## We assume that frames are in this order
    my @frames = (1, 2, 3, -1, -2, -3);
    my @peps = map( extract_pep($id, $lines[$_], $frames[$_/2]), (1,3,5,7,9,11) );
##    my @peps = map( extract_pep($id, $_), @lines[1,3,5,7,9,11] );
    ## this assumes that the frames are ordered as :
    ## +1, +2, +3, -1, -2, -3
    my @fwd = (1, @{$peps[0]});
    my @rev = (-1, @{$peps[3]});
    for my $i(1..2){
        @fwd = ($i+1, @{$peps[$i]}) if( length($peps[$i][0]) > length($fwd[1]) );
        @rev = (-1 * ($i+1), @{$peps[$i + 3]}) if( length($peps[$i+3][0]) > length($rev[1]) );
    }
    return($id, [@fwd], [@rev] );
}
sub calc_coverage {
    my( $nuc_length, $pep_length, $frame ) = @_;
    $frame = abs($frame);
    $nuc_length -= ($frame - 1);
    $nuc_length -= ($nuc_length % 3);
    return( 3 * $pep_length / $nuc_length );
}

## should be 12 lines per identifier
sub read_pep {
    my($file) = @_;
    my %pep = ();
    open(my $in, $file) || die "unable to open $file $!\n";
    while(<$in>){
        chomp;
        last if($_ !~ /^>\S+/);
        my @lines = ($_);
        for my $i(0..10){
            chomp($_ = <$in>);
            push @lines, $_;
        }
        my @orfs = extract_orfs(@lines);
        if(!defined($nuc{$orfs[0]})){
            die "no nucleotide sequence for : $orfs[0]\n";
        }
        @cov = map( { calc_coverage( length($nuc{$orfs[0]}), length( $$_[1] ), $$_[0]) } @orfs[1,2] );
        print "$orfs[0]\t", join("\t", @{$orfs[1]}), "\t$cov[0]\t", join("\t", @{$orfs[2]}), "\t$cov[1]\n";
        $pep{fwd}{length($orfs[1][1])}{$orfs[0]} = [(@{$orfs[1]}, $cov[0])];
        $pep{rev}{length($orfs[2][1])}{$orfs[0]} = [(@{$orfs[2]}, $cov[0])];
    }
    return(%pep);
}

sub read_fa {
    my($file) = @_;
    open(my $in, $file) || die "unable to open $file $!\n";
    my $id = "";
    my %seq = ();
    while(<$in>){
        chomp;
        if( $_ =~ /^>(\S+)\s*(.*)/ ){
            $id = $1;
            $seq{$id} = "";
            next;
        }
        $seq{$id} .= $_;
    }
    return( %seq );
}

sub rev_complement {
    my $seq = shift;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    $seq = reverse($seq);
    return($seq);
}

sub print_sequences {
    my($prefix, $hash_ref, $group_size, $min_length) = @_;
    my @lengths = sort {$b <=> $a} keys %{$hash_ref};
    my $of_count = 0;
    open(my $out, ">", sprintf("${prefix}_%04d.fa", $of_count / $group_size)) || die "unable to open sequence file for writin\
g: $prefix $of_count\n";
    for my $length(@lengths){
        last if $length < $min_length ;
        for my $id( keys %{$hash_ref->{$length}} ){
            $of_count++;
            if($of_count % $group_size == 0){
                open( $out, ">", sprintf("${prefix}_%04d.fa", $of_count / $group_size)) || die "unable to open sequence file \
for writing: $prefix $of_count\n";
            }
            print $out ">", $id, " $length ", join(":", @{$hash_ref->{$length}{$id}}[0,2,3,4,5]), "\n", $hash_ref->{$length}{\
$id}[1], "\n";
        }
    }
}

