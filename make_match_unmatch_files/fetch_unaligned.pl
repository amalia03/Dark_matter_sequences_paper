#!/bin/perl

($q_ids, $fq )=@ARGV;
$sam_name="$q_ids";
#print "Compiling indices for $sam_name\n";

open($fh, "<", $q_ids) || die " Could not open file $sam";

while(<$fh>){
    chomp;
    if($_=~ /^(\S+)/){
        $query{$1} = 1;
#       print $1, "\n";
    }
}

#print "Finished compiling indices\n";

open($fh2, "<", $fq) || die " Could not open file $fq\n";

$line_pos = 0;
$mapped = 0;
while(<$fh2>){
    $line_pos++;
    if ($line_pos % 2 == 1){
        if($_ =~ /^>(\S+)/){
            $mapped = defined( $query{$1} );
        }else{
            die "What the hell $!\n$_";
        }
    }
    print if ! $mapped;
}
print "\n";
