#!/bin/perl -w 

($fastafile, $sample_id) = @ARGV;

open($fh, "<", $fastafile) || die "Could not open file $fastafile/n $!";

while (<$fh>) {
    chomp;
    if($_ =~ /^>(.+)/ ){
        $id = $1;
	next;
    }    
    $nucl.=$_;
}
#print "This is the string: $nucl\n";

##If you place if instead of while, it will stop at the first instance, while takes in all values

while ($nucl=~/A/g) {$count_A++}
while ($nucl=~/G/g) {$count_G++}
while ($nucl=~/C/g) {$count_C++}
while ($nucl=~/T/g) {$count_T++}
$nucl_sum=$count_A+$count_T+$count_C+$count_G;
$perc_A=$count_A / $nucl_sum; 
$perc_T=$count_T / $nucl_sum; 
$perc_C=$count_C / $nucl_sum; 
$perc_G=$count_G / $nucl_sum; 

print sprintf("%.2f",$perc_A),"\t",sprintf("%.2f",$perc_C),"\t",sprintf("%.2f",$perc_T),"\t",sprintf("%.2f",$perc_G),"\n"; 

#print $nucl_sum, "\n";
#print "A nucleotides constitute ", sprintf("%.2f",$perc_A), " of the total nucleotides.\n"; 
#print "C nucleotides constitute ", sprintf("%.2f",$perc_C), " of the total nucleotides.\n"; 
#print "T nucleotides constitute ", sprintf("%.2f",$perc_T), " of the total nucleotides.\n"; 
#print "G nucleotides constitute ", sprintf("%.2f",$perc_G), " of the total nucleotides.\n"; 
