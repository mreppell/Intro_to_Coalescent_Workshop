#! usr/bin/perl 

use strict;
use Getopt::Long;

my $FULLFILE = "none";
my $opt = 0;
my %data;
my %sd_data;
my $OUT = "none";
my $output = 0;

GetOptions("f=s" => \$FULLFILE, "o=s" => \$OUT, "p=f"=>\$output);

if (($FULLFILE eq "none") or ($OUT eq "none")){
    die "Usage -f [Input Filename] -o [Output Filename] -p [Output format [0/1]]\n";
}

my @SHORTFILE = split /\//, $FULLFILE;
my $FILE = $SHORTFILE[$#SHORTFILE];

my $opt = 0;
foreach my $nval (split /\./, $FILE) {
    if ($nval eq "gz") {++$opt;}
}

if ($opt==0) {
    open(IN,$FULLFILE) or die "(1) Can't open $FULLFILE\n";
} else {
    open(IN,"gunzip -c $FULLFILE |") or die "(2) Can't open $FULLFILE\n";
}

open(OUT,">$OUT") or die "Can't open $OUT\n";

my $head1 = <IN>;
chomp($head1);
my @info1 =  split /\s+/, $head1;
my $prog = $info1[0];
my $haps = $info1[1];
my $runs = $info1[2]; 

for my $run (1..$runs) {
    for my $haplos (1..$haps) {
        $sd_data{$run}{$haplos} = 0;
        if ($run==1) {
            $data{$_} = 0;
        } 
    }
}

my $count = 0;
    GO: while (my $line = <IN>) {
        chomp($line);
        my @vals = split /\s+/, $line;
  #      print "$line\n";
        if ($vals[0] eq "segsites:") {
            $count++;
            my @sumarray;
            my $length = $vals[1];
            my $slength = $length-1;
            <IN>;
            my @overall;
            for (1..$length) {
                push(@overall, 0);
            }
            
            while (my $inline = <IN>) {
                    chomp($inline);
                    my @plow = split /\s+/, $inline;
                    
                    if (exists $plow[0]) {
                    my @in = split //, $plow[0];
                        
                        for (0..$slength) {
               #            print "$in[$_]\n";
                            $overall[$_] = $overall[$_] + $in[$_];
                        }
                    }
                    else {
                        for (1..$haps) {
                            my $spotcount = 0;
                            foreach my $spots (@overall) {
                                if ($spots == $_) {
                                    $spotcount++;
                                }
                            }
                         #   print "var: $_ $data{$_} + $spotcount\n";
                            $data{$_} = $data{$_} + $spotcount;
                            $sd_data{$count}{$_} = $spotcount;
                        }
                        next GO;
                    }
              }
             # print "here\n";
                   #     print "@overall\n";
                        for (1..$haps) {
                            my $spotcount = 0;
                            foreach my $spots (@overall) {
                                if ($spots == $_) {
                                    $spotcount++;
                                }
                            }
                        #    print "2 var: $_ $data{$_} + $spotcount\n";
                            $data{$_} = $data{$_} + $spotcount;
                            $sd_data{$count}{$_} = $spotcount;
                        }
                 
            
            
    }
}

#print "Count: $count\n";

if ($output==0) {
print OUT "ALLELE_COUNT AVG SD\n"; 
foreach my $num (sort {$a <=> $b} keys %data) {
    my $avg = $data{$num}/$runs;
    my $var = 0;
    foreach my $r (keys %sd_data) {
        $var = $var + ($sd_data{$r}{$num}-$avg)*($sd_data{$r}{$num}-$avg);
    }
    my $sd = sprintf("%.3f",sqrt((1/($runs-1))*$var));
    if ($num>0) {
        print OUT "$num $avg $sd\n";
    }
}
} 

if ($output==1) {
    foreach my $c (sort {$a <=> $b} keys %sd_data) {
        foreach my $v (1..$haps) {
            if (exists $sd_data{$c}{$v}) {
                print OUT "$sd_data{$c}{$v} ";
            } else {
                print OUT "0 ";
            }
        }
        print OUT "\n";
    }
}

#for my $run (keys %sd_data) {
#    for my $haplos (keys %{$sd_data{$run}}) {
#        print "Run$run $haplos $sd_data{$run}{$haplos}\n";
#    }
#} 
 
close(OUT);            
close(IN);
            
