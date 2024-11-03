#!/usr/bin/perl -w
#Import fasta file
($fasta_file, $fasta_l) = @ARGV;

#Read in the file
open($fh, "<", $fasta_file) || die "Could not open file $fasta_file/n $!";

while (<$fh>) {
    chomp;
    if($_ =~ /^>(\S+)/ ){
        $id = $1;
        $seqs{$id} = "";
        next;
    }
    $seqs{$id} .= $_;

}

$l_seq =length($seqs{$id});
#print "Length of fasta sequence : $l_seq","\n";

open($fh2, "<", $fasta_l) || die "Could not open file $fasta_file/n $!";
while(<$fh2>){
    chomp;
    @tmp = split /\t/, $_;
    $range = $tmp[1];
    $count++;
    $minimum= int(rand(length($seqs{$id})-$range));
#       print "Random minimum is : $minimum\n";
    $fa_sbstr  = substr $seqs{$id}, $minimum, $range;
    print ">ID$count\n$fa_sbstr\n";
}
