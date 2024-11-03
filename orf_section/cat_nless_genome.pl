#!/bin/perl -w

($fastafile, $group_name) = @ARGV;

open($fh, "<", $fastafile) || die "Could not open file $fastafile/n $!";

while (<$fh>) {
    chomp;
    if($_ =~ /^>(.+)/ ){
#        $id = $1;
        next;
    }
    $nucl .= $_ =~ s/N//gr;

}

print ">concatenated_",$group_name,"_genome","\n";

print "$nucl\n";
